clc
%format short
global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b

% Variabili simboliche:

% Ingressi Tau=[tau_phi tau_psi]
Tau=sym('tau',[2 1],'real');

% Disturbi di misura V=[v_psi v_phi_dot v_dx v_db];
V = sym('v',[4 1],'real');

% Disturbi di processo W = [w_tau_phi w_tau_psi]
W = sym('w',[2 1],'real');

% Vettore di stato X = [x y theta phi_dot psi psi_dot]
X = sym('x',[6 1],'real');

old_f = [X;Tau;V;W];
old_h = [X;V];

% Caricamento del dataset:
load('dataset')
dt = log_vars.dt;
t_max = log_vars.t_max;
meas_sens_psi = log_vars.meas_sens_psi;
meas_sens_phi_dot = log_vars.meas_sens_phi_dot;
meas_sens_dx = log_vars.meas_sens_dx;
meas_sens_db = log_vars.meas_sens_db;

x_hat = log_vars.X_hat;     % stima iniziale
P = log_vars.P;             % matrice di covarianza associata

Q = diag([log_vars.std_dev_tau_phi^2, log_vars.std_dev_tau_psi^2]); %matrice di disturbo di processo


% Funzioni simboliche per il calcolo della funzione di stato, del modello
% di misura e dei rispettivi jacobiani:
Jsym = getJacobian();

% Funzione di stato f simbolica
f = Jsym.f;
matlabFunction(f,'File','state_function_f_EKF','Vars',{old_f});
% Modello di misura h simbolico
h = Jsym.h;
matlabFunction(h,'File','observation_model_h_EKF','Vars',{old_h});
% Jacobiano simbolico di f rispetto allo stato X
F = Jsym.F;
matlabFunction(F,'File','jacobian_F','Vars',{old_f},'Optimize',false);
% Jacobiano simbolico di f rispetto al disturbo W
T = Jsym.T;
matlabFunction(T,'File','jacobian_T','Vars',{old_f});
% Jacobiano simbolico di h rispetto allo stato X
H = Jsym.H;
matlabFunction(H,'File','jacobian_H','Vars',{old_h});
% Jacobiano simbolico di h rispetto al disturbo V
M = Jsym.M;
matlabFunction(M,'File','jacobian_M','Vars',{old_h});


selection_vector = [false false false false]';  % seleziona quali misure sono state usate all'iterazione corrente
flag = [0 0 0 0]';  % tiene traccia dell'indice delle misure più recenti già utilizzate per ogni sensore
actual_meas = [0 0 0 0]';   % contiene i valori delle misure utilizzate all'iterazione corrente

% in = 1;   % secondo me possiamo toglierla questa variabile perché è
% esattamente uguale al k, quindi basta usare quello. Dentro getActualMeas
% gli sto passando k invece che in
k = 1;
for t = dt:dt:t_max
    %prediction step
    [x_hat(:,k+1), P] = prediction_EKF(x_hat(:,k), P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k));
    
    % restituisce ad ogni passo il vettore con le misure dei sensori
    % effettive e disponibili che non erano già state prese
    [actual_meas, selection_vector, flag] = getActualMeas(meas_sens_psi, meas_sens_phi_dot, meas_sens_dx, ...
                                                         meas_sens_db, flag, selection_vector, k, dt);

    % in = in +1;

    R = diag([log_vars.std_dev_psi^2, log_vars.std_dev_phi_dot^2, ...
          log_vars.std_dev_dx^2, log_vars.std_dev_db^2]); %matrice di errore di misura

    % correction step
    [x_hat(:,k+1), P, innovation(:,k)] = correction_EKF(x_hat(:,k+1), P, R, actual_meas, selection_vector);

    k = k + 1;
end

% Salvataggio nel dataset delle variabili per i plot
[x_estimation]=[x_hat(1,:)]';
[y_estimation]=[x_hat(2,:)]';
[theta_estimation] = [x_hat(3,:)]';
[phi_dot_estimation] = [x_hat(4,:)]';
[psi_estimation] = [x_hat(5,:)]';
[psi_dot_estimation] = [x_hat(6,:)]';

log_vars.x_estimation_EKF = x_estimation;
log_vars.y_estimation_EKF = y_estimation;
log_vars.theta_estimation_EKF = theta_estimation;
log_vars.phi_dot_estimation_EKF = phi_dot_estimation;
log_vars.psi_estimation_EKF = psi_estimation;
log_vars.psi_dot_estimation_EKF = psi_dot_estimation;
log_vars.innovation = innovation;

save('dataset','log_vars');

% Funzioni:

function  [x_hat, P] = prediction_EKF(x_hat, P, Q, dt, tau_phi, tau_psi)
    % Calcolo dei jacobiani numerici:
    new = [x_hat; tau_phi; tau_psi; 0; 0; 0; 0; 0; 0];
    % Jacobiani numerici calcolati con matlabFunction
    F = jacobian_F(new);
    T = jacobian_T(new);
%     % Jacobiani numerici calcolati con subs
%     F = double(subs(Jsym.F, old_f,new));
%     T = double(subs(Jsym.T, old_f, new));
    
    % Calcolo della funzione di stato al passo k
    % Funzione di stato numerica calcolata con matlabFunction
    f = state_function_f_EKF(new);
%     % Funzione di stato numerica calcolata con subs
%     f = double(subs(Jsym.f, old_f, new));
    
    % stima di X al passo k+1 noto k
    x_hat = x_hat + dt*f;
    x_hat(3) = wrapTo2Pi(x_hat(3));
    x_hat(5) = wrapToPi(x_hat(5));
    % Matrice di covarianza associata
    P = F*P*F' + T*Q*T';
end

function [x_hat, P, innovation] = correction_EKF(x_hat, P, R, actual_meas, selection_vector)
    % Calcolo dei jacobiani numerici:
    new = [x_hat; 0; 0; 0; 0];
    % Jacobiani numerici calcolati con matlabFunction
    H = jacobian_H(new);
    M = jacobian_M(new);
%     % Jacobiani numerici calcolati con subs
%     H = double(subs(Jsym.H, old_h, new));
%     M = double(subs(Jsym.M, old_h, new));

    % Calcolo del modello di osservazione al passo k
    % Modello di misura numerico calcolato con matlabFunction
    h = observation_model_h_EKF(new);
%     % Modello di misura numerico calcolato con subs
%     h = double(subs(Jsym.h, old_h, new));

    % Calcolo dell'innovazione
    innovation = [0 0 0 0]';
    counter = 0;
    if selection_vector(1) == false     % ridimensiona le matrici se la misura di psi non viene utilizzata
        h(1) = [];
        H(1,:) = [];
        M(1,:) = [];
        M(:,1) = [];
        R(1,:) = [];
        R(:,1) = [];
        innovation(1) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
        counter = counter + 1;
    end

    if selection_vector(2) == false     % ridimensiona le matrici se la misura di phi_dot non viene utilizzata 
        h(2-counter) = [];
        H(2-counter,:) = [];
        M(2-counter,:) = [];
        M(:,2-counter) = [];
        R(2-counter,:) = [];
        R(:,2-counter) = [];
        innovation(2) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
        counter = counter + 1;
    end
    
    if selection_vector(3) == false     % ridimensiona le matrici se la misura di dx non viene utilizzata
        h(3-counter) = [];
        H(3-counter,:) = [];
        M(3-counter,:) = [];
        M(:,3-counter) = [];
        R(3-counter,:) = [];
        R(:,3-counter) = [];
        innovation(3) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
        counter = counter + 1;
    end
    
    if selection_vector(4) == false     % ridimensiona le matrici se la misura di db non viene utilizzata
        h(4-counter) = [];
        H(4-counter,:) = [];
        M(4-counter,:) = [];
        M(:,4-counter) = [];
        R(4-counter,:) = [];
        R(:,4-counter) = [];
        innovation(4) = 100;    % 100 indica sul vettore innovation che la misura non è stata utilizzata
    end
    
    % innovazione attuale con solo le misure correnti fornite dai sensori
    e = actual_meas' - h;
    
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

    
    % Calcolo della covarianza associata all'innovazione
    S = H*P*H' + M*R*M';
    % Calcolo del guadagno di correzione
    L = P*H'*inv(S);
    
    % Stima statica BLUE
    if(isempty(e) == true)  % se non ci sono misure disponibili non viene fatta correzione
        x_hat = x_hat;
        P = P;
    else
        x_hat = x_hat + L*e;
        x_hat(3) = wrapTo2Pi(x_hat(3));
        x_hat(5) = wrapToPi(x_hat(5));
        % Calcolo della P con la forma di Joseph
        P = (eye(6) - L*H)*P*(eye(6) - L*H)' + L*M*R*M'*L';
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
