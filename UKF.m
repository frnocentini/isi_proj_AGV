clc
close all
global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b


%Tau=[tau_phi tau_psi]
Tau=sym('tau',[2 1],'real');

%V=[v_psi v_phi_dot v_dx v_db];
V = sym('v',[4 1],'real');

%W = [w_tau_phi w_tau_psi]
W = sym('w',[2 1],'real');

X = sym('x',[6 1],'real');

old_f = [X;W;Tau];
old_h = [X;V];

load('dataset')

x_hat = log_vars.X_hat;
P = log_vars.P;

Q = diag([log_vars.std_dev_tau_phi^2, log_vars.std_dev_tau_psi^2]); %matrice di disturbo di processo
R = diag([log_vars.std_dev_psi^2, log_vars.std_dev_phi_dot^2, ...
          log_vars.std_dev_dx^2, log_vars.std_dev_db^2]); %matrice di errore di misura

dt = log_vars.dt;
t_max = log_vars.t_max;
meas_sens_psi = log_vars.meas_sens_psi;
meas_sens_phi_dot = log_vars.meas_sens_phi_dot;
meas_sens_dx = log_vars.meas_sens_dx;
meas_sens_db = log_vars.meas_sens_db;

k = 1;
Jsym = getJacobian();
%f = Jsym.f;
%function_f = matlabFunction(f, 'File','state_function_f','Vars',{old_f});
%h = Jsym.h;
%function_h = matlabFunction(h);
for t = dt:dt:t_max

    %prediction step
    [x_hat(:,k+1), P] = prediction_UKF(x_hat(:,k), P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k), Jsym, old_f);
    %correction step
    [x_hat(:,k+1), P] = correction_UKF(x_hat(:,k+1), P, R, meas_sens_psi(k+1), meas_sens_phi_dot(k+1), meas_sens_dx(k+1), ...
        meas_sens_db(k+1), Jsym, old_h);

    k = k + 1;
end

% Salvataggio sul dataset delle variabili per i plot
[x_estimation]=[x_hat(1,:)]';
[y_estimation]=[x_hat(2,:)]';
[theta_estimation] = [x_hat(3,:)]';
[phi_dot_estimation] = [x_hat(4,:)]';
[psi_estimation] = [x_hat(5,:)]';
[psi_dot_estimation] = [x_hat(6,:)]';

log_vars.x_estimation_UKF = [x_estimation];
log_vars.y_estimation_UKF = [y_estimation];
log_vars.theta_estimation_UKF = [theta_estimation];
log_vars.phi_dot_estimation_UKF = [phi_dot_estimation];
log_vars.psi_estimation_UKF = [psi_estimation];
log_vars.psi_dot_estimation_UKF = [psi_dot_estimation];

save('dataset','log_vars');


function  [x_hat, P] = prediction_UKF(x_hat, P, Q, dt, tau_phi, tau_psi, Jsym, old_f)
    [x_hat, P] = UT_F([x_hat; 0; 0], blkdiag(P, Q), dt, tau_phi, tau_psi, Jsym, old_f);
end

function [sample_mean, sample_cov] = UT_F(prior_mean, prior_cov, dt, tau_phi, tau_psi, Jsym, old_f)
    
    % Parametri dell'algoritmo SUT che descrivono quanto sono scalati i
    % sigma point simmetrici rispetto al valore centrale
    alpha = 1;            % UT non scalata (lambda = k)
    beta = 2;             % ottimo nel caso gaussiano
    kappa = 0;            % serve per approssimare ordini superiori al secondo
    n = size(prior_mean,1); %dimensione del vettore media a priori
    lambda = alpha^2*(n+kappa)-n; %indica quanto distanti dal valore centrale bisogna prendere
                                  %i sigma points
    
    %Calcolo del vettore pesi per momento del primo ordine w_mean (vettore
    %riga)
    w0 = lambda/(n+lambda); % peso del sigma point centrale
    wi = 1/(2*(n+lambda)); %pesi dei sigma point simmetrici (con i != 0)
    w_mean(1) = w0;
    w_mean(2:2*n+1) = wi; 

    %Calcolo del vettore pesi per momenti del secondo ordine (w0_cov è w0'
    %della teoria) mentre gli altri si calcolano allo stesso modo
    w0_cov = w0+1-alpha^2+beta; % peso del sigma point centrale
    wi_cov = 1/(2*(n+lambda)); %pesi dei sigma point simmetrici (con i != 0)
    w_cov(1) = w0_cov;
    w_cov(2:2*n+1) = wi_cov;

    %Fattorizzazione della matrice di covarianza
    [U,S,~] = svd(prior_cov);
    GAMMA = U*S^(1/2);

    %Generazione dei sigma points
    sigma_points = prior_mean; %Il sigma point centrale corrisponde alla media a priori
    for i = 1:size(GAMMA,2)
        ai = sqrt(n/(1-w0));
        sigma_points(:,i+1)      = prior_mean + ai*GAMMA(:,i); %sigma point a destra di quello centrale (i>0)
        sigma_points(:,i+1+n) = prior_mean - ai*GAMMA(:,i); %sigma point a sinistra di quello centrale (i<0)
    end

    %Calcolo dei sigma points propagati
    for i = 1:size(sigma_points,2)
        new = [sigma_points(:,i); tau_phi; tau_psi];
        propagated_sigmapoints(1:6,i) = sigma_points(1:6,i) + dt*double(subs(Jsym.f, old_f, new)); 
    end

    %Media campionaria
    sample_mean = propagated_sigmapoints*w_mean';
    sample_mean(3,1) = wrapTo2Pi(atan2(sin(propagated_sigmapoints(3,:)*w_mean'),cos(propagated_sigmapoints(3,:)*w_mean')));
    sample_mean(5,1) = wrapToPi(atan2(sin(propagated_sigmapoints(5,:)*w_mean'),cos(propagated_sigmapoints(5,:)*w_mean')));

    tilde = propagated_sigmapoints - sample_mean;

    %Varianza campionaria
    sample_cov = zeros(6,6);
    for i = 1:size(sigma_points,2)
        sample_cov = sample_cov + w_cov(i)*tilde(:,i)*tilde(:,i)';
    end
end


function [x_hat, P] = correction_UKF(x_hat, P, R, meas_psi, meas_phi_dot, meas_dx, meas_db, Jsym, old_h)
    [hat_y, S, Pxy] = UT_H(x_hat, P, meas_psi, meas_phi_dot, meas_dx, meas_db, Jsym, old_h);
    
    % calcolo dell'innovazione
    y = [meas_psi, meas_phi_dot, meas_dx, meas_db]';
    e = y - hat_y;
    e(1) = wrapToPi(e(1,1))
    % calcolo della covarianza associata all'innovazione aggiungendo il
    % rumore
    S = S + R;
    % calcolo del guadagno di correzione
    L = Pxy*(S)^(-1);
    
    %stima statica BLUE
    x_hat = x_hat + L*e;
    x_hat(3) = wrapTo2Pi(x_hat(3));
    x_hat(5) = wrapToPi(x_hat(5));
    %calcolo della P con la forma di Joseph
    P = P - L*S*L';
    
end

function [pred_meas, pred_meas_cov, cross_cov] = UT_H(prior_mean, prior_cov, ...
    meas_psi, meas_phi_dot, meas_dx, meas_db, Jsym, old_h)
    % Parametri dell'algoritmo SUT che descrivono quanto sono scalati i
    % sigma point simmetrici rispetto al valore centrale
    alpha = 1;            % UT non scalata (lambda = k)
    beta = 2;             % ottimo nel caso gaussiano
    kappa = 0;            % serve per approssimare ordini superiori al secondo
    n = size(prior_mean,1); %dimensione del vettore media a priori
    lambda = alpha^2*(n+kappa)-n; %indica quanto distanti dal valore centrale bisogna prendere
                                  %i sigma points
    
    %Calcolo del vettore pesi per momento del primo ordine
    w0 = lambda/(n+lambda); % peso del sigma point centrale
    wi = 1/(2*(n+lambda)); %pesi dei sigma point simmetrici (con i != 0)
    w_mean(1) = w0;
    w_mean(2:2*n+1) = wi; 

    %Calcolo del vettore pesi per momenti del secondo ordine
    w0_cov = w0+1-alpha^2+beta; % peso del sigma point centrale
    wi_cov = 1/(2*(n+lambda)); %pesi dei sigma point simmetrici (con i != 0)
    w_cov(1) = w0_cov;
    w_cov(2:2*n+1) = wi_cov;

    %Fattorizzazione della matrice di covarianza
    [U,S,~] = svd(prior_cov);
    GAMMA = U*S^(1/2);

    %Generazione dei sigma points
    sigma_points = prior_mean; %Il sigma point centrale corrisponde alla media a priori
    for i = 1:size(GAMMA,2)
        ai = sqrt(n/(1-w0));
        sigma_points(:,i+1) = prior_mean + ai*GAMMA(:,i); %sigma point a destra di quello centrale (i>0)
        sigma_points(:,i+1+n) = prior_mean - ai*GAMMA(:,i); %sigma point a sinistra di quello centrale (i<0)
    end

    %Calcolo della funzione h sui sigma points
    for i = 1:size(sigma_points,2)
        new = [sigma_points(:,i); 0; 0; 0; 0];
        propagated_sigmapoints(:,i) = double(subs(Jsym.h, old_h, new));
    end

    %Media campionaria
    pred_meas = propagated_sigmapoints*w_mean';
    pred_meas(1,1) = wrapToPi(atan2(sin(propagated_sigmapoints(1,:)*w_mean'), cos(propagated_sigmapoints(1,:)*w_mean')));
    
    tilde = propagated_sigmapoints - pred_meas;
    tilde(1,:) = wrapToPi(atan2(sin(propagated_sigmapoints(1,:) - pred_meas(1,1)), cos(propagated_sigmapoints(1,:) - pred_meas(1,1))));

    %Varianza campionaria e covarianza incrociata
    pred_meas_cov = zeros(4,4);
    cross_cov = zeros(6,4);
    for i = 1:size(sigma_points,2)
        pred_meas_cov = pred_meas_cov + w_cov(i)*tilde(:,i)*tilde(:,i)';
        cross_cov = cross_cov + w_cov(i)*(sigma_points(:,i) - prior_mean)*tilde(:,i)';
    end
end