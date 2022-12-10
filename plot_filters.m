%% Caricamento del dataset
load('dataset')
x_real = log_vars.x_real;
y_real = log_vars.y_real;
theta_real = log_vars.theta_real;
phi_dot_real = log_vars.phi_dot_real;
psi_real = log_vars.psi_real;
psi_dot_real = log_vars.psi_dot_real;
dt = log_vars.dt;
%% Plot EKF

x_estimation_EKF = log_vars.x_estimation_EKF;
y_estimation_EKF = log_vars.y_estimation_EKF;
theta_estimation_EKF = log_vars.theta_estimation_EKF;
phi_dot_estimation_EKF = log_vars.phi_dot_estimation_EKF;
psi_estimation_EKF = log_vars.psi_estimation_EKF;
psi_dot_estimation_EKF = log_vars.psi_dot_estimation_EKF;


figure(1);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(x_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(x_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,400]); ylabel('x');
legend({'x real','x stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(y_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(y_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,400]); ylabel('y');
legend({'y reale','y stimata'},'orientation','horizontal','location','southoutside');
title(t1,'EKF results'); xlabel(t1,'Time (sec)');


figure(2);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(theta_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(theta_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,400]); ylabel('theta');
legend({'theta reale','theta stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(phi_dot_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(phi_dot_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,400]); ylabel('phi_dot');
legend({'phi\_dot reale','phi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t2,'EKF results'); xlabel(t2,'Time (sec)');


figure(3);
t3 = tiledlayout(2,1);
nexttile; 
% Plot psi
plot(psi_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(psi_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,400]); ylabel('psi');
legend({'psi reale','psi stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot psi_dot
plot(psi_dot_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(psi_dot_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,400]); ylabel('psi_dot');
legend({'psi\_dot reale','psi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t3,'EKF results'); xlabel(t3,'Time (sec)');



%% Plot UKF

x_estimation_UKF = log_vars.x_estimation_UKF;
y_estimation_UKF = log_vars.y_estimation_UKF;
theta_estimation_UKF = log_vars.theta_estimation_UKF;
phi_dot_estimation_UKF = log_vars.phi_dot_estimation_UKF;
psi_estimation_UKF = log_vars.psi_estimation_UKF;
psi_dot_estimation_UKF = log_vars.psi_dot_estimation_UKF;


figure(1);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(x_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on; 
plot(x_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]); ylabel('x'); 
legend({'x reale','x stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(y_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(y_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]); ylabel('y');
legend({'y reale','y stimata'},'orientation','horizontal','location','southoutside');
title(t1,'UKF results'); xlabel(t1,'Time (sec)');


figure(2);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(theta_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(theta_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]); ylabel('theta');
legend({'theta reale','theta stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(phi_dot_real,'LineWidth',1,color = '[0.9290 0.6940 0.1250]'); hold on;
plot(phi_dot_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]); ylabel('phi_dot');
legend({'phi\_dot reale','phi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t2,'UKF results'); xlabel(t2,'Time (sec)');


figure(3);
t3 = tiledlayout(2,1);
nexttile; 
% Plot psi
plot(psi_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(psi_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]); ylabel('psi');
legend({'psi reale','psi stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot psi_dot
plot(psi_dot_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on; 
plot(psi_dot_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]); ylabel('psi_dot');
legend({'psi\_dot reale','psi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t3,'UKF results'); xlabel(t3,'Time (sec)');

%% Plot della traiettoria risultante dallo stimatore EKF
f = 3;
figure(1)
hold on
plot(x_estimation_EKF', y_estimation_EKF')

for k = 1 : size(x_estimation_EKF, 1)
    x = x_estimation_EKF(k,:);
    th = theta_estimation_EKF(k,:);
    y = y_estimation_EKF(k,:);
    h = plot([x, x + f*cos(th)],[y, y + f*sin(th)], 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    
    pause(dt)
    delete(h)
end


