%%% Rauch–Tung–Striebel Smoother
clc

load('dataset')
load('EKF_struct.mat')

x_real = log_vars.x_real;
y_real = log_vars.y_real;

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
    x_hat_smooth_theta(k)= log_EKF(:,k).x_hat_smoothed(3);
    x_hat_smooth_psi(k)= log_EKF(:,k).x_hat_smoothed(5);
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

%% Animazione
x_hat_smooth_x = x_hat_smooth_x';
x_hat_smooth_y = x_hat_smooth_y';
x_hat_smooth_theta = x_hat_smooth_theta';
x_hat_smooth_psi = x_hat_smooth_psi';

close all
clc
figure(1)
axis equal
hold on
plot(x_hat_smooth_x', x_hat_smooth_y')
plot(x_real', y_real', '--');
coord_BEx = [d,-L-0.5];
coord_BEy = [d+0.1,d+0.1];
coord_CDx = [d,-L-0.5];
coord_CDy = [-d-0.1,-d-0.1];
coord_BCx = [d,d];
coord_BCy = [d+0.1,-d-0.1];
coord_DEx = [-L-0.5,-L-0.5];
coord_DEy = [-d-0.1,+d+0.1];
coord_wax_1 = [-rp,rp];
coord_way_1 = [d,d];
coord_wax_2 = [rp, -rp];
coord_way_2 = [-d,-d];

scale_factor = 5;
xlim([min(x_hat_smooth_x-10) max(x_hat_smooth_x+10)])
ylim([min(x_hat_smooth_y-10) max(x_hat_smooth_y+10)])

%animazione dell'AGV lungo la traiettoria e della ruota posteriore
for k = 1 : size(x_hat_smooth_x, 1)
    x = x_hat_smooth_x(k,:);
    th = x_hat_smooth_theta(k,:);
    y = x_hat_smooth_y(k,:);
    psi = x_hat_smooth_psi(k,:);
    R_zt = [cos(th) , -sin(th); sin(th), cos(th)];
    BE_R = scale_factor*R_zt*[coord_BEx; coord_BEy];
    BC_R = scale_factor*R_zt*[coord_BCx; coord_BCy];
    CD_R = scale_factor*R_zt*[coord_CDx; coord_CDy];
    DE_R = scale_factor*R_zt*[coord_DEx; coord_DEy];
    wp = scale_factor*R_zt*[-L+rp*cos(psi), -L-rp*cos(psi);rp*sin(psi), -rp*sin(psi)];
    wa1 = scale_factor*R_zt*[coord_wax_1; coord_way_1];
    wa2 = scale_factor*R_zt*[coord_wax_2; coord_way_2];
        
    BE = plot(x+BE_R(1,:),y+BE_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    BC = plot(x+BC_R(1,:),y+BC_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    CD= plot(x+CD_R(1,:),y+CD_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    DE = plot(x+DE_R(1,:),y+DE_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');   
    w_p = plot(x+wp(1,:),y+wp(2,:),'LineWidth', 2, color = 'r');
    w_a1 = plot(x+wa1(1,:),y+wa1(2,:),'LineWidth', 2, color = 'r');
    w_a2 = plot(x+wa2(1,:),y+wa2(2,:),'LineWidth', 2, color = 'r');

    title('AGV plot')
    legend('regolarizzata', 'groundtruth')
    pause(dt)
    delete(BE)
    delete(BC)
    delete(CD)
    delete(DE)
    delete(w_p)
    delete(w_a1)
    delete(w_a2)
end

hold off

