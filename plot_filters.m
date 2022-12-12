%% Caricamento del dataset
global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b

load('dataset')
x_real = log_vars.x_real;
y_real = log_vars.y_real;
theta_real = log_vars.theta_real;
phi_dot_real = log_vars.phi_dot_real;
psi_real = log_vars.psi_real;
psi_dot_real = log_vars.psi_dot_real;
innovation_EKF = log_vars.innovation;
%innovation_UKF = log_vars.innovation_UKF;

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

%plot innovazione lungo tutte le componenti

figure(4);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(innovation_EKF(1,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
xlim([0,400]); 
legend('Innovazione per psi','orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(innovation_EKF(2,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on;
xlim([0,400]);
legend('Innovazione per phi\_dot','orientation','horizontal','location','southoutside');
title(t1,'Innovazione EKF'); xlabel(t1,'Time (sec)');


figure(5);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(innovation_EKF(3,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
xlim([0,400]); 
legend({'Innovazione per dx'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(innovation_EKF(4,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
xlim([0,400]); 
legend('Innovazione per db','southoutside');
title(t2,'Innovazione EKF'); xlabel(t2,'Time (sec)');

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

%plot innovazione lungo tutte le componenti

figure(4);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(innovation_UKF(1,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
xlim([0,400]); 
legend('Innovazione per psi','orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(innovation_UKF(2,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on;
xlim([0,400]);
legend('Innovazione per phi\_dot','orientation','horizontal','location','southoutside');
title(t1,'Innovazione EKF'); xlabel(t1,'Time (sec)');


figure(5);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(innovation_UKF(3,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
xlim([0,400]); 
legend({'Innovazione per dx'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(innovation_UKF(4,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
xlim([0,400]); 
legend('Innovazione per db','southoutside');
title(t2,'Innovazione EKF'); xlabel(t2,'Time (sec)');

%% Plot della traiettoria risultante dallo stimatore EKF
close all
figure(1)
hold on
plot(x_estimation_EKF', y_estimation_EKF')
plot(x_real', y_real', '--');
coord_BEx = [d,-L-0.5];
coord_BEy = [d+0.1,d+0.1];
coord_CDx = [d,-L-0.5];
coord_CDy = [-d-0.1,-d-0.1];
coord_BCx = [d,d];
coord_BCy = [d+0.1,-d-0.1];
coord_DEx = [-L-0.5,-L-0.5];
coord_DEy = [-d-0.1,+d+0.1];

xlim([min(x_estimation_EKF-10) max(x_estimation_EKF+10)])
ylim([min(y_estimation_EKF-10) max(y_estimation_EKF+10)])

%animazione dell'AGV lungo la traiettoria e della ruota posteriore
for k = 1 : size(x_estimation_EKF, 1)
    x = x_estimation_EKF(k,:);
    th = theta_estimation_EKF(k,:);
    y = y_estimation_EKF(k,:);
    psi = psi_estimation_EKF(k,:);
    R_zt = [cos(th) , -sin(th); sin(th), cos(th)];
    BE_R = R_zt*[coord_BEx; coord_BEy];
    BC_R = R_zt*[coord_BCx; coord_BCy];
    CD_R = R_zt*[coord_CDx; coord_CDy];
    DE_R = R_zt*[coord_DEx; coord_DEy];
    wp = R_zt*[-L+rp*2*cos(psi), -L-rp*2*cos(psi);rp*2*sin(psi), -rp*2*sin(psi)];
    BE = plot(x+BE_R(1,:),y+BE_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    BC = plot(x+BC_R(1,:),y+BC_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    CD= plot(x+CD_R(1,:),y+CD_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    DE = plot(x+DE_R(1,:),y+DE_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');   
    w_p = plot(x+wp(1,:),y+wp(2,:),'LineWidth', 2, color = 'r');
    title('AGV plot')
    legend('stimata', 'groundtruth')
    pause(dt)
    delete(BE)
    delete(BC)
    delete(CD)
    delete(DE)
    delete(w_p)
end


%% Plot della traiettoria risultante dallo stimatore UKF
close all
figure(1)
hold on
axis equal
plot(x_estimation_UKF', y_estimation_UKF')
plot(x_real', y_real', '--');
coord_BEx = [d,-L-0.5];
coord_BEy = [d+0.1,d+0.1];
coord_CDx = [d,-L-0.5];
coord_CDy = [-d-0.1,-d-0.1];
coord_BCx = [d,d];
coord_BCy = [d+0.1,-d-0.1];
coord_DEx = [-L-0.5,-L-0.5];
coord_DEy = [-d-0.1,+d+0.1];

xlim([min(x_estimation_UKF-10) max(x_estimation_UKF+10)])
ylim([min(y_estimation_UKF-10) max(y_estimation_UKF+10)])

%animazione dell'AGV lungo la traiettoria e della ruota posteriore
for k = 1 : size(x_estimation_UKF, 1)
    x = x_estimation_UKF(k,:);
    th = theta_estimation_UKF(k,:);
    y = y_estimation_UKF(k,:);
    psi = psi_estimation_UKF(k,:);
    R_zt = [cos(th) , -sin(th); sin(th), cos(th)];
    BE_R = R_zt*[coord_BEx; coord_BEy];
    BC_R = R_zt*[coord_BCx; coord_BCy];
    CD_R = R_zt*[coord_CDx; coord_CDy];
    DE_R = R_zt*[coord_DEx; coord_DEy];
    wp = R_zt*[-L+rp*2*cos(psi), -L-rp*2*cos(psi);rp*2*sin(psi), -rp*2*sin(psi)];
    BE = plot(x+BE_R(1,:),y+BE_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    BC = plot(x+BC_R(1,:),y+BC_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    CD= plot(x+CD_R(1,:),y+CD_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');
    DE = plot(x+DE_R(1,:),y+DE_R(2,:), 'LineWidth', 2, color = '[0.4940 0.1840 0.5560]');   
    w_p = plot(x+wp(1,:),y+wp(2,:),'LineWidth', 2, color = 'r');
    title('AGV plot')
    legend('stimata', 'groundtruth')
    pause(dt)
    delete(BE)
    delete(BC)
    delete(CD)
    delete(DE)
    delete(w_p)
end
