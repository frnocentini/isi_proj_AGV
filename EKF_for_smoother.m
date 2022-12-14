clc

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
R = diag([log_vars.std_dev_psi^2, log_vars.std_dev_phi_dot^2, ...
          log_vars.std_dev_dx^2, log_vars.std_dev_db^2]); %matrice di errore di misura

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


k = 0;
for t = dt:dt:t_max
    
    k = k+1;
    %prediction step
    [x_hat, P, F] = prediction_EKF(x_hat, P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k));
    log_EKF(k).x_hat_pred= x_hat;
    log_EKF(k).P_pred = P;
    log_EKF(k).F_matrix = F;

    %correction step
    [x_hat, P] = correction_EKF(x_hat, P, R, meas_sens_psi(k+1), meas_sens_phi_dot(k+1), meas_sens_dx(k+1), ...
                                        meas_sens_db(k+1));

    log_EKF(k).x_hat_corr= x_hat;
    log_EKF(k).P_corr = P;

    
end

  save('EKF_struct', 'log_EKF');
% Salvataggio nel dataset delle variabili per i plot
% [x_estimation]=[x_hat(1,:)]';
% [y_estimation]=[x_hat(2,:)]';
% [theta_estimation] = [x_hat(3,:)]';
% [phi_dot_estimation] = [x_hat(4,:)]';
% [psi_estimation] = [x_hat(5,:)]';
% [psi_dot_estimation] = [x_hat(6,:)]';
% 
% log_vars.x_estimation_EKF = x_estimation;
% log_vars.y_estimation_EKF = y_estimation;
% log_vars.theta_estimation_EKF = theta_estimation;
% log_vars.phi_dot_estimation_EKF = phi_dot_estimation;
% log_vars.psi_estimation_EKF = psi_estimation;
% log_vars.psi_dot_estimation_EKF = psi_dot_estimation;
% log_vars.innovation = [e];

% save('dataset','log_vars');

% Funzioni:
function  [x_hat, P, F] = prediction_EKF(x_hat, P, Q, dt, tau_phi, tau_psi)
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

function [x_hat, P] = correction_EKF(x_hat, P, R, meas_psi, meas_phi_dot, meas_dx, meas_db)
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
    h = observation_model_h_UKF(new);
%     % Modello di misura numerico calcolato con subs
%     h = double(subs(Jsym.h, old_h, new));

    % Calcolo dell'innovazione
    e = [meas_psi; meas_phi_dot; meas_dx; meas_db] - h;
    e(1) = wrapToPi(e(1,1))
    
    % Calcolo della covarianza associata all'innovazione
    S = H*P*H' + M*R*M';
    % Calcolo del guadagno di correzione
    L = P*H'*inv(S);

    % Stima statica BLUE
    x_hat = x_hat + L*e;
    x_hat(3) = wrapTo2Pi(x_hat(3));
    x_hat(5) = wrapToPi(x_hat(5));
    % Calcolo della P con la forma di Joseph
    P = (eye(6) - L*H)*P*(eye(6) - L*H)' + L*M*R*M'*L';
end