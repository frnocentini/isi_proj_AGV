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

x_hat = log_vars.X_hat;     % stima iniziale
P = log_vars.P;             % matrice di covarianza associata

Q = diag([log_vars.std_dev_tau_phi^2, log_vars.std_dev_tau_psi^2]); %matrice di disturbo di processo
R = diag([log_vars.std_dev_psi^2, log_vars.std_dev_phi_dot^2, ...
          log_vars.std_dev_dx^2, log_vars.std_dev_db^2]); %matrice di errore di misura

dt = log_vars.dt;
t_max = log_vars.t_max;
meas_sens_psi = log_vars.meas_sens_psi;
meas_sens_phi_dot = log_vars.meas_sens_phi_dot;
meas_sens_dx = log_vars.meas_sens_dx;
meas_sens_db = log_vars.meas_sens_db;
% x_real = log_vars.x;
% y_real = log_vars.y;
% theta_real = log_vars.theta;
% phi_dot_real = log_vars.phi_dot;
% psi_real = log_vars.psi;
% psi_dot_real = log_vars.psi_dot;

% Funzioni simboliche per il calcolo della funzione di stato, del modello
% di misura e dei rispettivi jacobiani:
Jsym = getJacobian();

% % Funzione di stato f simbolica
% f = Jsym.f;
% matlabFunction(f,'File','state_function_f','Vars',{old_f});
% % Modello di misura h simbolico
% h = Jsym.h;
% matlabFunction(h,'File','observation_model_h','Vars',{old_h});
% % Jacobiano simbolico di f rispetto allo stato X
% F = Jsym.F;
% matlabFunction(F,'File','function_F','Vars',{old_f});
% % Jacobiano simbolico di f rispetto al disturbo W
% T = Jsym.T;
% matlabFunction(T,'File','function_T','Vars',{old_f});
% % Jacobiano simbolico di h rispetto allo stato X
% H = Jsym.H;
% matlabFunction(H,'File','function_H','Vars',{old_h});
% % Jacobiano simbolico di h rispetto al disturbo V
% M = Jsym.M;
% matlabFunction(M,'File','function_M','Vars',{old_h});


k = 1;
for t = dt:dt:t_max

    %prediction step
    [x_hat(:,k+1), P] = prediction_EKF(x_hat(:,k), P, Q, dt, log_vars.tau_phi(k), log_vars.tau_psi(k), Jsym, old_f);
%     [fun(k)]=f(1);

    %correction step
    [x_hat(:,k+1), P,e,h] = correction_EKF(x_hat(:,k+1), P, R, meas_sens_psi(k), meas_sens_phi_dot(k), meas_sens_dx(k), meas_sens_db(k), Jsym, old_h);
%     [x_stimato(k)]=x_hat(1);
%     [y_stimato(k)]=x_hat(2);
%     [theta_stimato(k)]=x_hat(3);
%     [phi_dot_stimato(k)] = x_hat(4);
%     [psi_stimato(k)] = x_hat(5);
%     [psi_dot_stimato(k)] = x_hat(6);
%     [oss(k)]=h(3);

    k = k + 1;
end

% [x_stima]=[x_stimato]';
% [y_stima]=[y_stimato]';
% [theta_stima] = [theta_stimato]';
% [phi_dot_stima] = [phi_dot_stimato]';
% [psi_stima] = [psi_stimato]';
% [psi_dot_stima] = [psi_dot_stimato]';
% 
% subplot(6,1,1)
% plot(x_real,'b');
% hold on;
% plot(x_stima,'r');
% hold off;
% legend('xreale','xstimata');
% xlim([0,400]); ylim([0,20]);
% 
% subplot(6,1,2)
% plot(y_real,'b');
% hold on;
% plot(y_stima,'r');
% hold off;
% legend('yreale','ystimata');
% xlim([0,400]); ylim([0,20]);
% 
% subplot(6,1,3)
% plot(theta_real,'b');
% hold on;
% plot(theta_stima,'r');
% hold off;
% legend('thetaReale','thetaStimato');
% xlim([0,400]); ylim([0,20]);
% 
% subplot(6,1,4)
% plot(phi_dot_real,'b');
% hold on;
% plot(phi_dot_stima,'r');
% hold off;
% legend('phi_dotreale','phi_dotstimata');
% xlim([0,400]); ylim([0,20]);
% 
% subplot(6,1,5)
% plot(psi_real,'b');
% hold on;
% plot(psi_stima,'r');
% hold off;
% legend('psireale','psistimata');
% xlim([0,400]); ylim([0,20]);
% 
% subplot(6,1,6)
% plot(psi_dot_real,'b');
% hold on;
% plot(psi_dot_stima,'r');
% hold off;
% legend('psidotReale','psidotStimata');
% xlim([0,400]); ylim([0,20]);

function  [x_hat, P] = prediction_EKF(x_hat, P, Q, dt, tau_phi, tau_psi, Jsym, old_f)
    % Calcolo dei jacobiani numerici:
    new = [x_hat; tau_phi; tau_psi; 0; 0; 0; 0; 0; 0];
    %F = function_F(new);
    F = double(subs(Jsym.F, old_f,new));
    %T = function_T(new);
    T = double(subs(Jsym.T, old_f, new));
    
    % Calcolo della funzione di stato al passo k
    %f = state_function_f(new);
    f = double(subs(Jsym.f, old_f, new));
    
    % stima di X al passo k+1 noto k
    x_hat = x_hat + dt*f;
    x_hat(3) = wrapTo2Pi(x_hat(3));
    x_hat(5) = wrapToPi(x_hat(5));
    % Matrice di covarianza associata
    P = F*P*F' + T*Q*T';
end

function [x_hat, P, e,h] = correction_EKF(x_hat, P, R, meas_psi, meas_phi_dot, meas_dx, meas_db, Jsym, old_h)
    % Calcolo dei jacobiani numerici:
    new = [x_hat; 0; 0; 0; 0];
    %H = function_H(new);
    H = double(subs(Jsym.H, old_h, new));
    %M = function_M(new);
    M = double(subs(Jsym.M, old_h, new));

    % Calcolo del modello di osservazione al passo k
    %h = observation_model_h(new);
    h = double(subs(Jsym.h, old_h, new));

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
