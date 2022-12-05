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
J = getJacobian();

for t = dt:dt:t_max

    k = k + 1;
    %prediction step
    [x_hat, P] = prediction_EKF(x_hat, P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k), J, old_f);
    %correction step
    [x_hat, P] = correction_EKF(x_hat, P, R, meas_sens_psi(k), meas_sens_phi_dot(k), meas_sens_dx(k), meas_sens_db(k), J, old_h);
    
end

function  [x_hat, P] = prediction_EKF(x_hat, P, Q, dt, tau_phi, tau_psi, J, old_f)
    new = [x_hat; tau_phi; tau_psi; 0; 0; 0; 0; 0; 0]';

    F = double(subs(J.F, old_f, new));
    T = double(subs(J.T, old_f, new));
    
    f = double(subs(J.f, old_f, new));
    
    x_hat = x_hat + dt*f;
    P = F*P*F' + T*Q*T';
    %D = 
    %f = X + log_vars.dt*f;

end

function [x_hat, P] = correction_EKF(x_hat, P, R, meas_psi, meas_phi_dot, meas_dx, meas_db, J, old_h)
    new = [x_hat; 0; 0; 0; 0]';
    H = double(subs(J.H, old_h, new));
    M = double(subs(J.M, old_h, new));
    h = double(subs(J.h, old_h, new));
    
    %calcolo dell'innovazione
    e = [meas_psi; meas_phi_dot; meas_dx; meas_db] - h;
    %calcolo della covarianza associata all'innovazione
    S = H*P*H' + M*R*M';

    %calcolo del guadagno di correzione
    L = P*H'*inv(S)
    %stima statica BLUE
    x_hat = x_hat + L*e;
    %calcolo della P con la forma di Joseph
    P = (eye(6) - L*H)*P*(eye(6) - L*H)' + L*M*R*M'*L';
    
end