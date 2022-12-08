clc

global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b
syms x y theta phi_dot psi_ psi_dot tau_phi_s tau_psi_s ...
       w_psi_ w_phi_dot w_dx w_db w_tau_psi w_tau_phi

old_f = [x y theta phi_dot psi_ psi_dot tau_phi_s tau_psi_s w_psi_ w_phi_dot w_dx w_db w_tau_psi w_tau_phi];
old_h = [x y theta phi_dot psi_ psi_dot w_psi_ w_phi_dot w_dx w_db];


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

k = 0;
Jsym = getJacobian();
f = Jsym.f;
function_f = matlabFunction(f);
h = Jsym.h;
function_h = matlabFunction(h);
for t = dt:dt:t_max

    k = k + 1;
    %prediction step
    [x_hat, P] = prediction_UKF(x_hat, P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k), function_f);
    %correction step
    [x_hat, P,e] = correction_UKF(x_hat, P, R, meas_sens_psi(k), meas_sens_phi_dot(k), meas_sens_dx(k), ...
        meas_sens_db(k), function_h);
%     [x_stimato(k)]=x_hat(1);
%     [oss(k)]=h(3);
end

function  [x_hat, P] = prediction_UKF(x_hat, P, Q, dt, tau_phi, tau_psi, function_f)
    [x_hat, P] = UT_F([x_hat; 0; 0], blkdiag(P, Q), dt, tau_phi, tau_psi, function_f);
end

function [sample_mean, sample_cov] = UT_F(prior_mean, prior_cov, dt, tau_phi, tau_psi, function_f)

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
        sigma_points(:,i+1)      = prior_mean + ai*GAMMA(:,i); %sigma point a destra di quello centrale (i>0)
        sigma_points(:,i+1+enne) = prior_mean - ai*GAMMA(:,i); %sigma point a sinistra di quello centrale (i<0)
    end

    %Calcolo dei sigma points propagati
    for i = 1:size(sigma_points,2)
        propagated_sigmapoints(:,i) = function_f(sigma_points(4,i),sigma_points(5,i),...
            sigma_points(6,i),tau_phi,tau_psi,sigma_points(3,i),sigma_points(7,i),sigma_points(8,i)); 
    end

    %Media campionaria
    sample_mean = propagated_sigmapoints*w_mean';
    sample_mean(2,1) = atan2(sin(propagated_sigmapoints(2,:)*w_mean'),cos(propagated_sigmapoints(2,:)*w_mean'));
    sample_mean(6,1) = atan2(sin(propagated_sigmapoints(6,:)*w_mean'),cos(propagated_sigmapoints(6,:)*w_mean'));

    tilde = propagated_sigmapoints - sample_mean;

    %Varianza campionaria
    sample_cov = zeros(6,6);
    for i = 1:size(sigma_points,2)
        sample_cov = sample_cov + w_cov(i)*tilde(:,i)*tilde(:,i)';
    end
end


function [x_hat, P, e] = correction_UKF(x_hat, P, R, meas_psi, meas_phi_dot, meas_dx, meas_db, function_h)
    [hat_y, S, Pxy] = UT_H(x_hat, P, meas_psi, meas_phi_dot, meas_dx, meas_db,function_h);
    
    % calcolo dell'innovazione
    y = [meas_psi, meas_phi_dot, meas_dx, meas_db];
    e = y - hat_y;
    e(1) = wrapToPi(e(1,1));
    % calcolo della covarianza associata all'innovazione aggiungendo il
    % rumore
    S = S + R;
    % calcolo del guadagno di correzione
    L = Pxy*inv(S);
    
    %stima statica BLUE
    x_hat = x_hat + L*e;
    %calcolo della P con la forma di Joseph
    P = P - L*S*L';
    
end

function [pred_meas, pred_meas_cov, cross_cov] = UT_H(prior_mean, prior_cov, ...
    meas_psi, meas_phi_dot, meas_dx, meas_db, function_h)
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
        sigma_points(:,i+1)      = prior_mean + ai*GAMMA(:,i); %sigma point a destra di quello centrale (i>0)
        sigma_points(:,i+1+enne) = prior_mean - ai*GAMMA(:,i); %sigma point a sinistra di quello centrale (i<0)
    end

    %Calcolo della funzione h sui sigma points
    for i = 1:size(sigma_points,2)
        propagated_sigmapoints(:,i) = function_h(sigma_points(4,i),sigma_points(5,i),...
            0,0,0,0,sigma_points(1,i),sigma_points(2,i));
    end

    %Media campionaria
    pred_meas = propagated_sigmapoints*w_mean';
    pred_meas(2,1) = atan2(sin(propagated_sigmapoints(2,:)*w_mean'),cos(propagated_sigmapoints(2,:)*w_mean'));
    
    tilde = propagated_sigmapoints - pred_meas;
    tilde(2,1) = atan2(sin(propagated_sigmapoints - pred_meas),cos(propagated_sigmapoints - pred_meas));

    %Varianza campionaria e covarianza incrociata
    pred_meas_cov = zeros(6,6);
    for i = 1:size(sigma_points,2)
        pred_meas_cov = pred_meas_cov + w_cov(i)*tilde(:,i)*tilde(:,i)';
        cross_cov = cross_cov + w_cov(i)*(sigma_points(:,i) - pred_meas)*tilde(:,i)';
    end
end