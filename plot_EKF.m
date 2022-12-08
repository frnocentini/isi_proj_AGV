%% Plot

% Caricamento del dataset
load('dataset')

x_real = log_vars.x_real;
y_real = log_vars.y_real;
theta_real = log_vars.theta_real;
phi_dot_real = log_vars.phi_dot_real;
psi_real = log_vars.psi_real;
psi_dot_real = log_vars.psi_dot_real;

x_estimation = log_vars.x_estimation;
y_estimation = log_vars.y_estimation;
theta_estimation = log_vars.theta_estimation;
phi_dot_estimation = log_vars.phi_dot_estimation;
psi_estimation = log_vars.psi_estimation;
psi_dot_estimation = log_vars.psi_dot_estimation;


figure(1);
t1 = tiledlayout(2,1);
nexttile; plot(x_real,'b'); hold on; plot(x_estimation,'r');
xlim([0,400]);
ylabel('x');
legend({'xreale','xstimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(y_real,'b'); hold on; plot(y_estimation,'r');
xlim([0,400]);
ylabel('y');
title(t1,'EKF results'); xlabel(t1,'Time (sec)');
legend({'yreale','ystimata'},'orientation','horizontal','location','southoutside');

figure(2);
t2 = tiledlayout(2,1);
nexttile; plot(theta_real,'b'); hold on; plot(theta_estimation,'r');
xlim([0,400]);
ylabel('theta');
legend({'thetareale','thetastimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(phi_dot_real,'b'); hold on; plot(phi_dot_estimation,'r');
xlim([0,400]);
ylabel('phi_dot');
title(t2,'EKF results'); xlabel(t2,'Time (sec)');
legend({'phi_dotreale','phi_dotstimata'},'orientation','horizontal','location','southoutside');

figure(3);
t3 = tiledlayout(2,1);
nexttile; plot(psi_real,'b'); hold on; plot(psi_estimation,'r');
xlim([0,400]);
ylabel('psi');
legend({'psireale','psistimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(psi_dot_real,'b'); hold on; plot(psi_dot_estimation,'r');
xlim([0,400]);
ylabel('psi_dot');
title(t3,'EKF results'); xlabel(t3,'Time (sec)');
legend({'psi_dotreale','psi_dotstimata'},'orientation','horizontal','location','southoutside');
