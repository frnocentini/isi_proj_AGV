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

dt = log_vars.dt;
t_max = log_vars.t_max;
meas_sens_psi = log_vars.meas_sens_psi;
meas_sens_phi_dot = log_vars.meas_sens_phi_dot;
meas_sens_dx = log_vars.meas_sens_dx;
meas_sens_db = log_vars.meas_sens_db;

% Funzioni simboliche per il calcolo della funzione di stato e del modello
% di misura:
Jsym = getJacobian();
% Funzione di stato f simbolica
f = Jsym.f;
matlabFunction(f,'File','state_function_f_UKF','Vars',{old_f});
% Modello di misura h simbolico
h = Jsym.h;
matlabFunction(h,'File','observation_model_h_UKF','Vars',{old_h});

selection_vector = [false false false false]';  % seleziona quali misure sono state usate all'iterazione corrente
flag = [0 0 0 0]';  % tiene traccia dell'indice delle misure più recenti già utilizzate per ogni sensore
actual_meas = [0 0 0 0]';   % contiene i valori delle misure utilizzate all'iterazione corrente

k = 1;
for t = dt:dt:t_max

    %prediction step
    [x_hat(:,k+1), P] = prediction_UKF(x_hat(:,k), P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k)); 

    % restituisce ad ogni passo il vettore con le misure dei sensori
    % effettive e disponibili che non erano già state prese
    [actual_meas, selection_vector, flag] = getActualMeas(meas_sens_psi, meas_sens_phi_dot, meas_sens_dx, ...
                                                         meas_sens_db, flag, selection_vector, k, dt);

    R = diag([log_vars.std_dev_psi^2, log_vars.std_dev_phi_dot^2, ...
          log_vars.std_dev_dx^2, log_vars.std_dev_db^2]);   % matrice di errore di misura

    %correction step
    [x_hat(:,k+1), P, innovation(:,k+1)] = correction_UKF(x_hat(:,k+1), P, R, actual_meas, selection_vector);

    k = k + 1;
end

% Salvataggio sul dataset delle variabili per i plot
[x_estimation]=[x_hat(1,:)]';
[y_estimation]=[x_hat(2,:)]';
[theta_estimation] = [x_hat(3,:)]';
[phi_dot_estimation] = [x_hat(4,:)]';
[psi_estimation] = [x_hat(5,:)]';
[psi_dot_estimation] = [x_hat(6,:)]';

log_vars.x_estimation_UKF = x_estimation;
log_vars.y_estimation_UKF = y_estimation;
log_vars.theta_estimation_UKF = theta_estimation;
log_vars.phi_dot_estimation_UKF = phi_dot_estimation;
log_vars.psi_estimation_UKF = psi_estimation;
log_vars.psi_dot_estimation_UKF = psi_dot_estimation;
log_vars.innovation_UKF = innovation;

save('dataset','log_vars');

% Funzioni:

function  [x_hat, P] = prediction_UKF(x_hat, P, Q, dt, tau_phi, tau_psi)
    [x_hat, P] = UT_F([x_hat; 0; 0], blkdiag(P, Q), dt, tau_phi, tau_psi);
end


function [sample_mean, sample_cov] = UT_F(prior_mean, prior_cov, dt, tau_phi, tau_psi)
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
        sigma_points(:,i+1) = prior_mean + ai*GAMMA(:,i); %sigma point a destra di quello centrale (i>0)
        sigma_points(:,i+1+n) = prior_mean - ai*GAMMA(:,i); %sigma point a sinistra di quello centrale (i<0)
    end

    %Calcolo dei sigma points propagati
    for i = 1:size(sigma_points,2)
        new = [sigma_points(:,i); tau_phi; tau_psi];
        %propagated_sigmapoints(1:6,i) = sigma_points(1:6,i) + dt*double(subs(Jsym.f, old_f, new)); % calcolo con subs
        propagated_sigmapoints(1:6,i) = sigma_points(1:6,i) + dt*state_function_f_UKF(new); % calcolo con matlabFunction
    end

    %Media campionaria
    sample_mean = propagated_sigmapoints*w_mean';
    sample_mean(3,1) = wrapTo2Pi(atan2(sin(propagated_sigmapoints(3,:))*w_mean',cos(propagated_sigmapoints(3,:))*w_mean'));
    sample_mean(5,1) = wrapToPi(atan2(sin(propagated_sigmapoints(5,:))*w_mean',cos(propagated_sigmapoints(5,:))*w_mean'));

    tilde = propagated_sigmapoints - sample_mean;

    %Varianza campionaria
    sample_cov = zeros(6,6);
    for i = 1:size(sigma_points,2)
        sample_cov = sample_cov + w_cov(i)*tilde(:,i)*tilde(:,i)';
    end
end


function [x_hat, P, innovation] = correction_UKF(x_hat, P, R, actual_meas, selection_vector)
    [hat_y, S, Pxy] = UT_H(x_hat, P);
    
    % calcolo dell'innovazione
    innovation = [0 0 0 0]';
    counter = 0;
    if selection_vector(1) == false     % ridimensiona le matrici se la misura di psi non viene utilizzata
        hat_y(1) = [];
        S(1,:) = [];
        S(:,1) = [];
        R(1,:) = [];
        R(:,1) = [];
        Pxy(:,1) = [];
        innovation(1) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
        counter = counter + 1;
    end

    if selection_vector(2) == false     % ridimensiona le matrici se la misura di phi_dot non viene utilizzata 
        hat_y(2-counter) = [];
        S(2-counter,:) = [];
        S(:,2-counter) = [];
        R(2-counter,:) = [];
        R(:,2-counter) = [];
        Pxy(:,2-counter) = [];
        innovation(2) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
        counter = counter + 1;
    end
    
    if selection_vector(3) == false     % ridimensiona le matrici se la misura di dx non viene utilizzata
        hat_y(3-counter) = [];
        S(3-counter,:) = [];
        S(:,3-counter) = [];
        R(3-counter,:) = [];
        R(:,3-counter) = [];
        Pxy(:,3-counter) = [];
        innovation(3) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
        counter = counter + 1;
    end
    
    if selection_vector(4) == false     % ridimensiona le matrici se la misura di db non viene utilizzata
        hat_y(4-counter) = [];
        S(4-counter,:) = [];
        S(:,4-counter) = [];
        R(4-counter,:) = [];
        R(:,4-counter) = [];
        Pxy(:,4-counter) = [];
        innovation(4) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
    end

    % innovazione attuale con solo le misure correnti fornite dai sensori
    e = actual_meas' - hat_y;
    
    counter = 0;
    if selection_vector(1) == true
        e(1) = wrapToPi(e(1,1));
        counter = counter + 1;
        innovation(1) = e(counter);
    end

    if selection_vector(2) == true
        counter = counter + 1;
        innovation(2) = e(counter);
    end

    if selection_vector(3) == true
        counter = counter + 1;
        innovation(3) = e(counter);
    end

    if selection_vector(4) == true
        counter = counter + 1;
        innovation(4) = e(counter);
    end

    % Calcolo della covarianza associata all'innovazione aggiungendo il rumore
    S = S + R;
    % Calcolo del guadagno di correzione
    L = Pxy*(S)^(-1);
    
    % Stima statica BLUE
    if(isempty(e) == true)  % se non ci sono misure disponibili non viene fatta correzione
        x_hat = x_hat;
        P = P;
    else
        x_hat = x_hat + L*e;
        x_hat(3) = wrapTo2Pi(x_hat(3));
        x_hat(5) = wrapToPi(x_hat(5));
        % Calcolo della P
        P = P - L*S*L';
    end    
end

function [pred_meas, pred_meas_cov, cross_cov] = UT_H(prior_mean, prior_cov)
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
        %propagated_sigmapoints(:,i) = double(subs(Jsym.h, old_h, new)); % calcolo con subs
        propagated_sigmapoints(:,i) = observation_model_h_UKF(new); % calcolo con matlabFunction
    end

    %Media campionaria
    pred_meas = propagated_sigmapoints*w_mean';
    pred_meas(1,1) = wrapToPi(atan2(sin(propagated_sigmapoints(1,:))*w_mean', cos(propagated_sigmapoints(1,:))*w_mean'));
    
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


function [meas, selection_vector, flag] = getActualMeas(meas_sens_psi, meas_sens_phi_dot, ...
                                            meas_sens_dx, meas_sens_db, flag, selection_vector, in, dt)
    meas = [];  % vettore dei valori delle misure attuali
    
    % all'interno dei cicli while si ricava l'ultima misura disponibile
    % fornita dai sensori; se non c'è nessuna misura allora non si aggiorna
    % la variabile flag

    % MISURA DI PSI:
    count = 0;  % indica se sono entrato nel while
    while(((flag(1)) < size(meas_sens_psi.data,1)) && (meas_sens_psi.time(flag(1)+1) <= dt*in))
        count = count + 1;
        flag(1) = flag(1) + 1;
    end
    
    count_size_meas = 0;
    if(count == 0)
        selection_vector(1) = false;    % non c'è nessuna misura disponibile
    else
        count_size_meas = count_size_meas + 1;
        selection_vector(1) = true;     % la misura è disponibile
        meas(count_size_meas) = meas_sens_psi.data(flag(1));    % salvo la misura su meas[]
    end
    
    % MISURA DI PHI_DOT
    count = 0;  % indica se sono entrato nel while
    while(((flag(2)) < size(meas_sens_phi_dot.data,1)) && (meas_sens_phi_dot.time(flag(2)+1) <= dt*in))
        count = count + 1;
        flag(2) = flag(2) + 1;
    end

    if(count == 0)
        selection_vector(2) = false;    % non c'è nessuna misura disponibile
    else
        count_size_meas = count_size_meas + 1;
        selection_vector(2) = true;     % la misura è disponibile
        meas(count_size_meas) = meas_sens_phi_dot.data(flag(2));    % salvo la misura su meas[]
    end
    
    % MISURA DI DX
    count = 0;  % indica se sono entrato nel while
    while(((flag(3)) < size(meas_sens_dx.data,1)) && (meas_sens_dx.time(flag(3)+1) <= dt*in))
        count = count + 1;
        flag(3) = flag(3) + 1;
    end

    if(count == 0)
        selection_vector(3) = false;    % non c'è nessuna misura disponibile
    else
        count_size_meas = count_size_meas + 1;
        selection_vector(3) = true;     % la misura è disponibile
        meas(count_size_meas) = meas_sens_dx.data(flag(3));     % salvo la misura su meas[]
    end
    
    % MISURA DI DB
    count = 0;  % indica se sono entrato nel while
    while(((flag(4)) < size(meas_sens_db.data,1)) && (meas_sens_db.time(flag(4)+1) <= dt*in))
        count = count + 1;
        flag(4) = flag(4) + 1;
    end


    if(count == 0)
        selection_vector(4) = false;    % non c'è nessuna misura disponibile
    else
        count_size_meas = count_size_meas + 1;
        selection_vector(4) = true;     % la misura è disponibile
        meas(count_size_meas) = meas_sens_db.data(flag(4));     % salvo la misura su meas[]
    end
end