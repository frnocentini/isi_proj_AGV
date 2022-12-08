% Caricamento del dataset
load('dataset')
x_real = log_vars.x_real;
y_real = log_vars.y_real;
theta_real = log_vars.theta_real;
phi_dot_real = log_vars.phi_dot_real;
psi_real = log_vars.psi_real;
psi_dot_real = log_vars.psi_dot_real;

%% Plot EKF

x_estimation_EKF = log_vars.x_estimation_EKF;
y_estimation_EKF = log_vars.y_estimation_EKF;
theta_estimation_EKF = log_vars.theta_estimation_EKF;
phi_dot_estimation_EKF = log_vars.phi_dot_estimation_EKF;
psi_estimation_EKF = log_vars.psi_estimation_EKF;
psi_dot_estimation_EKF = log_vars.psi_dot_estimation_EKF;


figure(1);
t1 = tiledlayout(2,1);
nexttile; plot(x_real,'b'); hold on; plot(x_estimation_EKF,'r');
xlim([0,400]);
ylabel('x');
legend({'xreale','xstimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(y_real,'b'); hold on; plot(y_estimation_EKF,'r');
xlim([0,400]);
ylabel('y');
title(t1,'EKF results'); xlabel(t1,'Time (sec)');
legend({'yreale','ystimata'},'orientation','horizontal','location','southoutside');

figure(2);
t2 = tiledlayout(2,1);
nexttile; plot(theta_real,'b'); hold on; plot(theta_estimation_EKF,'r');
xlim([0,400]);
ylabel('theta');
legend({'thetareale','thetastimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(phi_dot_real,'b'); hold on; plot(phi_dot_estimation_EKF,'r');
xlim([0,400]);
ylabel('phi_dot');
title(t2,'EKF results'); xlabel(t2,'Time (sec)');
legend({'phi_dotreale','phi_dot stimata'},'orientation','horizontal','location','southoutside');

figure(3);
t3 = tiledlayout(2,1);
nexttile; plot(psi_real,'b'); hold on; plot(psi_estimation_EKF,'r');
xlim([0,400]);
ylabel('psi');
legend({'psireale','psistimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(psi_dot_real,'b'); hold on; plot(psi_dot_estimation_EKF,'r');
xlim([0,400]);
ylabel('psi_dot');
title(t3,'EKF results'); xlabel(t3,'Time (sec)');
legend({'psi\_dotreale','psi_dotstimata'},'orientation','horizontal','location','southoutside');

%% Plot UKF

x_estimation_UKF = log_vars.x_estimation_UKF;
y_estimation_UKF = log_vars.y_estimation_UKF;
theta_estimation_UKF = log_vars.theta_estimation_UKF;
phi_dot_estimation_UKF = log_vars.phi_dot_estimation_UKF;
psi_estimation_UKF = log_vars.psi_estimation_UKF;
psi_dot_estimation_UKF = log_vars.psi_dot_estimation_UKF;


figure(1);
t1 = tiledlayout(2,1);
nexttile; plot(x_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on; 
plot(x_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]);
ylabel('x');
legend({'x reale','x stimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(y_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(y_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]);
ylabel('y');
title(t1,'UKF results'); xlabel(t1,'Time (sec)');
legend({'y reale','y stimata'},'orientation','horizontal','location','southoutside');

figure(2);
t2 = tiledlayout(2,1);
nexttile; plot(theta_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(theta_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]);
ylabel('theta');
legend({'theta reale','theta stimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(phi_dot_real,'LineWidth',1,color = '[0.9290 0.6940 0.1250]'); hold on;
plot(phi_dot_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]);
ylabel('phi_dot');
title(t2,'UKF results'); xlabel(t2,'Time (sec)');
legend({'phi\_dot reale','phi\_dot stimata'},'orientation','horizontal','location','southoutside');

figure(3);
t3 = tiledlayout(2,1);
nexttile; plot(psi_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(psi_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]);
ylabel('psi');
legend({'psi reale','psi stimata'},'orientation','horizontal','location','southoutside');
nexttile; plot(psi_dot_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on; 
plot(psi_dot_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,400]);
ylabel('psi_dot');
title(t3,'UKF results'); xlabel(t3,'Time (sec)');
legend({'psi\_dot reale','psi\_dot stimata'},'orientation','horizontal','location','southoutside');

