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


xlabel('x [m]')
ylabel('y [m]')

ground_truth = plot(log_vars.x_real', log_vars.y_real', 'g--', 'LineWidth', 1.2);
for k = 1:1:size(log_EKF,2)
         
    x_hat_corr_x(k)= log_EKF(:,k).x_hat_corr(1);
    x_hat_corr_y(k)= log_EKF(:,k).x_hat_corr(2);
    x_hat_smooth_x(k)= log_EKF(:,k).x_hat_smoothed(1);
    x_hat_smooth_y(k)= log_EKF(:,k).x_hat_smoothed(2);
%     x_hat_corr_tot_in = x_hat_corr_tot;
%     x_hat_smooth_tot = [log_EKF(:,k).x_hat_smoothed];
end


xlim([min(x_hat_smooth_x) max(x_hat_smooth_x)])
ylim([min(x_hat_smooth_y) max(x_hat_smooth_y)])

estimated_plot = plot(x_hat_corr_x, x_hat_corr_y,'r', 'LineWidth', 1.2);
smoothed_plot = plot(x_hat_smooth_x, x_hat_smooth_y,'b', 'LineWidth', 1.2);
% 
% 
% ground_init_truth = plot(log_vars.s_Tx(1), log_vars.s_Ty(1), 'go', 'MarkerFaceColor', 'g');
% estimated_init_plot = plot(log_EKF.x_hat_correction(1,1), log_EKF.x_hat_correction(1,2),'ro', 'MarkerFaceColor', 'r');
% smoothed_init_plot = plot(log_EKF.x_hat_smoothed(1,1), log_EKF.x_hat_smoothed(1,2),'bo', 'MarkerFaceColor', 'b');
% 
legend('GT','EKF','RTS')
