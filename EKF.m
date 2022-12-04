clc

global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b
syms x y theta phi_dot psi_ psi_dot tau_phi_s tau_psi_s ...
       w_psi_ w_phi_dot w_dx w_db w_tau_psi w_tau_phi

old = [x y theta phi_dot psi_ psi_dot tau_phi_s tau_psi_s w_psi_ w_phi_dot w_dx w_db w_tau_psi w_tau_phi];


load('dataset')

x_hat = log_vars.X_hat;
P = log_vars.P;

Q = diag([log_vars.std_dev_tau_phi^2, log_vars.std_dev_tau_psi^2]); %matrice di disturbo di processo
R = diag([log_vars.std_dev_psi^2, log_vars.std_dev_phi_dot^2, ...
          log_vars.std_dev_dx^2, log_vars.std_dev_db^2]); %matrice di errore di misura

dt = log_vars.dt;
t_max = log_vars.t_max;

k = 0;
J = getJacobian();

for t = dt:dt:t_max

    k = k + 1;
    %prediction step
    [x_hat, P] = prediction_EKF(x_hat, P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k), J, old);
    %correction step
    
end

function  [x_hat, P] = prediction_EKF(x_hat, P, Q, dt, tau_phi, tau_psi, J, old)
    new = [x_hat; tau_phi; tau_psi; 0; 0; 0; 0; 0; 0]';
    F = subs(J.F, old, new);
    T = subs(J.T, old, new);
    H = subs(J.H, old, new);
    M = subs(J.M, old, new);
    
    x_hat
    P = F*P*F' + T*Q*T';
    %D = 
    %f = X + log_vars.dt*f;

end