%%% Rauch–Tung–Striebel Smoother
clc

load('dataset')
load('EKF_struct.mat')

% La procedura in avanti viene svolta nel file EKF_smoother.mat dove si
% ricavano x_hat(k|k), x_hat(k+1|k), P(k|k), P(k+1|k)

%x_hat(n|n)
log_EKF(size(log_EKF,2)).x_hat_smoothed = log_EKF(size(log_EKF,2)).x_hat_corr;

%P(n|n)
log_EKF(size(log_EKF,2)).P_smoothed = log_EKF(size(log_EKF,2)).P_corr;

%Procedura in indietro
for k = size(log_EKF,2)-1:-1:1
    Ck = log_EKF(k).P_corr*(log_EKF(k+1).F_matrix)'*(log_EKF(k+1).P_pred)^-1;

    log_EKF(k).x_hat_smoothed = log_EKF(k).x_hat_corr + Ck* ... 
        (log_EKF(k+1).x_hat_smoothed - log_EKF(k+1).x_hat_pred);

    log_EKF(k).P_smoothed = log_EKF(k).P_corr + Ck* ... 
        (log_EKF(k+1).P_smoothed - log_EKF(k+1).P_pred)*Ck';


end

figure(152)
clf
hold on
grid on
box on
axis equal
xlim([-50 150])
ylim([-50 150])
xlabel('x [m]')
ylabel('y [m]')

ground_truth = plot(log_vars.x_real', log_vars.y_real', 'g--', 'LineWidth', 1.2);
x_hat_plot = log_EKF(.x_hat_corr;
estimated_plot = plot(log_EKF.x_hat_corr, log_EKF.x_hat_correction(:,2),'r', 'LineWidth', 1.2);
smoothed_plot = plot(log_EKF.x_hat_smoothed(:,1), log_EKF.x_hat_smoothed(:,2),'b', 'LineWidth', 1.2);


ground_init_truth = plot(log_vars.s_Tx(1), log_vars.s_Ty(1), 'go', 'MarkerFaceColor', 'g');
estimated_init_plot = plot(log_EKF.x_hat_correction(1,1), log_EKF.x_hat_correction(1,2),'ro', 'MarkerFaceColor', 'r');
smoothed_init_plot = plot(log_EKF.x_hat_smoothed(1,1), log_EKF.x_hat_smoothed(1,2),'bo', 'MarkerFaceColor', 'b');

legend('GT','EKF','RTS')
