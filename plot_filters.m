%% Caricamento del dataset
clc
global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b

load('dataset')
x_real = log_vars.x_real;
y_real = log_vars.y_real;
theta_real = log_vars.theta_real;
phi_dot_real = log_vars.phi_dot_real;
psi_real = log_vars.psi_real;
psi_dot_real = log_vars.psi_dot_real;
t_max = log_vars.t_max;
dt = log_vars.dt;
innovation_EKF = log_vars.innovation_EKF;
%innovation_UKF = log_vars.innovation_UKF;
dt = log_vars.dt;

x_estimation_EKF = log_vars.x_estimation_EKF;
y_estimation_EKF = log_vars.y_estimation_EKF;
theta_estimation_EKF = log_vars.theta_estimation_EKF;
phi_dot_estimation_EKF = log_vars.phi_dot_estimation_EKF;
psi_estimation_EKF = log_vars.psi_estimation_EKF;
psi_dot_estimation_EKF = log_vars.psi_dot_estimation_EKF;


% x_estimation_UKF = log_vars.x_estimation_UKF;
% y_estimation_UKF = log_vars.y_estimation_UKF;
% theta_estimation_UKF = log_vars.theta_estimation_UKF;
% phi_dot_estimation_UKF = log_vars.phi_dot_estimation_UKF;
% psi_estimation_UKF = log_vars.psi_estimation_UKF;
% psi_dot_estimation_UKF = log_vars.psi_dot_estimation_UKF;

%% Plot EKF

figure(1);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(x_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(x_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,t_max/dt]); ylabel('x');
legend({'x real','x stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(y_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(y_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,t_max/dt]); ylabel('y');
legend({'y reale','y stimata'},'orientation','horizontal','location','southoutside');
title(t1,'EKF results'); xlabel(t1, 'Iterations');


figure(2);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(theta_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(theta_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,t_max/dt]); ylabel('theta');
legend({'theta reale','theta stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(phi_dot_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(phi_dot_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,t_max/dt]); ylabel('phi\_dot');
legend({'phi\_dot reale','phi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t2,'EKF results'); xlabel(t2, 'Iterations');


figure(3);
t3 = tiledlayout(2,1);
nexttile; 
% Plot psi
plot(psi_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(psi_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,t_max/dt]); ylabel('psi');
legend({'psi reale','psi stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot psi_dot
plot(psi_dot_real,'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(psi_dot_estimation_EKF,'LineWidth',1, color = '[0.6350 0.0780 0.1840]');
xlim([0,t_max/dt]); ylabel('psi_dot');
legend({'psi\_dot reale','psi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t3,'EKF results'); xlabel(t3, 'Iterations');

% Plot innovazione lungo tutte le componenti
figure(4);
t1 = tiledlayout(2,1);
nexttile; 
% Plot innovazione per psi
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(1,i) == 100
        plot(t,innovation_EKF(1,i),' ');hold on;
    else
        plot(t,innovation_EKF(1,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
%plot(innovation_EKF(1,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
xlim([0,t_max]); 
legend('Innovazione per psi','orientation','horizontal','location','southoutside');

nexttile; 
% Plot innovazione per phi\_dot
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(2,i) == 100
        plot(t,innovation_EKF(2,i),' ');hold on;
    else
        plot(t,innovation_EKF(2,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); 
%plot(innovation_EKF(2,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on;
legend('Innovazione per phi\_dot','orientation','horizontal','location','southoutside');
title(t1,'Innovazione EKF'); xlabel(t1,'Time (s)');


figure(5);
t2 = tiledlayout(2,1);
nexttile; 
% Plot Innovazione per dx
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(3,i) == 100
        plot(t,innovation_EKF(3,i),' ');hold on;
    else
        plot(t,innovation_EKF(3,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); ylim([-1,1]);
%plot(innovation_EKF(3,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
legend({'Innovazione per dx'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot innovazione per db
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(4,i) == 100
        plot(t,innovation_EKF(4,i),' ');hold on;
    else
        plot(t,innovation_EKF(4,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); ylim([-1,1]);
%plot(innovation_EKF(4,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
legend('Innovazione per db','orientation','horizontal','location','southoutside');
title(t2,'Innovazione EKF'); xlabel(t2,'Time (s)');

% Plot dell'errore di stima (differenza tra stato stimato e stato reale)
estimation_error(1,:) = [x_real - x_estimation_EKF]';
estimation_error(2,:) = [y_real - y_estimation_EKF]';
estimation_error(3,:) = [theta_real - theta_estimation_EKF]';
estimation_error(4,:) = [phi_dot_real - phi_dot_estimation_EKF]';
estimation_error(5,:) = [psi_real - psi_estimation_EKF];
estimation_error(6,:) = [psi_dot_real - psi_dot_estimation_EKF];

figure(6);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(estimation_error(1,:),'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('estimation error for x');

nexttile; 
% Plot y
plot(estimation_error(2,:),'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('estimation error for y');
title(t1,'EKF estimation error'); xlabel(t1, 'Iterations');


figure(7);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(estimation_error(3,:),'LineWidth',1, color = '[0.4940 0.1840 0.5560]'); 
xlim([0,t_max/dt]); ylabel('estimation error for theta');

nexttile; 
% Plot phi_dot
plot(estimation_error(4,:),'LineWidth',1, color = '[0.4940 0.1840 0.5560]');  
xlim([0,t_max/dt]); ylabel('estimation error for phi_dot');
title(t2,'EKF estimation error'); xlabel(t2, 'Iterations');


figure(8);
t3 = tiledlayout(2,1);
nexttile; 
% Plot psi
plot(estimation_error(5,:),'LineWidth',1, color = '[0.4940 0.1840 0.5560]'); 
xlim([0,t_max/dt]); ylabel('estimation error for phi_dot');

nexttile; 
% Plot psi_dot
plot(estimation_error(6,:),'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('estimation error for phi_dot');
title(t3,'EKF estimation error'); xlabel(t3, 'Iterations');


%% Plot UKF


figure(1);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(x_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on; 
plot(x_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('x'); 
legend({'x reale','x stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(y_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(y_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('y');
legend({'y reale','y stimata'},'orientation','horizontal','location','southoutside');
title(t1,'UKF results'); xlabel(t1, 'Iterations');


figure(2);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(theta_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(theta_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('theta');
legend({'theta reale','theta stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(phi_dot_real,'LineWidth',1,color = '[0.9290 0.6940 0.1250]'); hold on;
plot(phi_dot_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('phi_dot');
legend({'phi\_dot reale','phi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t2,'UKF results'); xlabel(t2, 'Iterations');


figure(3);
t3 = tiledlayout(2,1);
nexttile; 
% Plot psi
plot(psi_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on;
plot(psi_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('psi');
legend({'psi reale','psi stimata'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot psi_dot
plot(psi_dot_real,'LineWidth',1, color = '[0.9290 0.6940 0.1250]'); hold on; 
plot(psi_dot_estimation_UKF,'LineWidth',1, color = '[0.4940 0.1840 0.5560]');
xlim([0,t_max/dt]); ylabel('psi_dot');
legend({'psi\_dot reale','psi\_dot stimata'},'orientation','horizontal','location','southoutside');
title(t3,'UKF results'); xlabel(t3, 'Iterations');

% Plot innovazione lungo tutte le componenti
figure(4);
t1 = tiledlayout(2,1);
nexttile; 
% Plot innovazione per psi
t = 0;
for i = 1 : size(innovation_UKF,2)
    t = t + dt;
    if innovation_UKF(1,i) == 100
        plot(t,innovation_UKF(1,i),' ');hold on;
    else
        plot(t,innovation_UKF(1,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); 
legend('Innovazione per psi','orientation','horizontal','location','southoutside');

nexttile; 
% Plot innovazione per phi\_dot
t = 0;
for i = 1 : size(innovation_UKF,2)
    t = t + dt;
    if innovation_UKF(2,i) == 100
        plot(t,innovation_UKF(2,i),' ');hold on;
    else
        plot(t,innovation_UKF(2,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); 
legend('Innovazione per phi\_dot','orientation','horizontal','location','southoutside');
title(t1,'Innovazione UKF'); xlabel(t1,'Time (s)');


figure(5);
t2 = tiledlayout(2,1);
nexttile; 
% Plot Innovazione per dx
t = 0;
for i = 1 : size(innovation_UKF,2)
    t = t + dt;
    if innovation_UKF(3,i) == 100
        plot(t,innovation_UKF(3,i),' ');hold on;
    else
        plot(t,innovation_UKF(3,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); ylim([-1,1]);
legend({'Innovazione per dx'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot innovazione per db
t = 0;
for i = 1 : size(innovation_UKF,2)
    t = t + dt;
    if innovation_UKF(4,i) == 100
        plot(t,innovation_UKF(4,i),' ');hold on;
    else
        plot(t,innovation_UKF(4,i),'o','LineWidth',1, color = '[0.3010 0.7450 0.9330]');hold on;
    end
end
xlim([0,t_max]); ylim([-1,1]);
legend('Innovazione per db','orientation','horizontal','location','southoutside');
title(t2,'Innovazione UKF'); xlabel(t2,'Time (s)');

%% Plot della traiettoria risultante dallo stimatore EKF
close all
clc
figure(1)
axis equal
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
coord_wax_1 = [-rp,rp];
coord_way_1 = [d,d];
coord_wax_2 = [rp, -rp];
coord_way_2 = [-d,-d];

scale_factor = 5;
xlim([min(x_estimation_EKF-10) max(x_estimation_EKF+10)])
ylim([min(y_estimation_EKF-10) max(y_estimation_EKF+10)])

%animazione dell'AGV lungo la traiettoria e della ruota posteriore
for k = 1 : size(x_estimation_EKF, 1)
    x = x_estimation_EKF(k,:);
    th = theta_estimation_EKF(k,:);
    y = y_estimation_EKF(k,:);
    psi = psi_estimation_EKF(k,:);
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
    legend('stimata', 'groundtruth')
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

%% confronto delle traiettorie
plot(x_estimation_EKF', y_estimation_EKF')
plot(x_estimation_UKF', y_estimation_UKF')
plot(x_real', y_real', '--');
legend('stima EKF','stima UKF', 'groundtruth');
title(t2,'Innovazione UKF'); xlabel(t2,'Time (s)');
%% Plot della traiettoria risultante dallo stimatore UKF
close all
figure(1)
hold on
axis equal
plot(x_estimation_UKF', y_estimation_UKF')

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

xlim([min(x_estimation_UKF-10) max(x_estimation_UKF+10)])
ylim([min(y_estimation_UKF-10) max(y_estimation_UKF+10)])

%animazione dell'AGV lungo la traiettoria e della ruota posteriore
for k = 1 : size(x_estimation_UKF, 1)
    x = x_estimation_UKF(k,:);
    th = theta_estimation_UKF(k,:);
    y = y_estimation_UKF(k,:);
    psi = psi_estimation_UKF(k,:);
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
    legend('stimata', 'groundtruth')
    pause(dt)
    delete(BE)
    delete(BC)
    delete(CD)
    delete(DE)
    delete(w_p)
    delete(w_a1)
    delete(w_a2)
end

%% confronto delle traiettorie
plot(x_estimation_EKF', y_estimation_EKF')
plot(x_estimation_UKF', y_estimation_UKF')
plot(x_real', y_real', '--');
legend('stima EKF','stima UKF', 'groundtruth');
title(t2,'Innovazione UKF'); xlabel(t2,'Time (s)');

%% confronto tra filtri componente per componente

figure(1);
t1 = tiledlayout(2,1);
nexttile; 
% Plot x
plot(x_real,'--','LineWidth',0.5, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(x_estimation_EKF,'LineWidth',0.5, color = '[0.6350 0.0780 0.1840]'); hold on;
plot(x_estimation_UKF,'LineWidth',0.5, color = 'g');
xlim([0,t_max/dt]); ylabel('x');
legend({'x real','x stimata EKF','x stimata UKF'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot y
plot(y_real,'--','LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(y_estimation_EKF,'LineWidth',0.5, color = '[0.6350 0.0780 0.1840]');hold on;
plot(y_estimation_UKF,'LineWidth',0.5, color = 'g');

xlim([0,t_max/dt]); ylabel('y');
legend({'y reale','y stimata EKF','y stimata UKF' },'orientation','horizontal','location','southoutside');
title(t1,'EKF results'); xlabel(t1, 'Iterations');


figure(2);
t2 = tiledlayout(2,1);
nexttile; 
% Plot theta
plot(theta_real,'--','LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(theta_estimation_EKF,'LineWidth',0.5, color = '[0.6350 0.0780 0.1840]');hold on;
plot(theta_estimation_UKF,'LineWidth',0.5, color = 'g'); 
xlim([0,t_max/dt]); ylabel('theta');
legend({'theta reale','theta stimata EKF', 'theta stimata UKF'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot phi_dot
plot(phi_dot_real,'--','LineWidth',0.5, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(phi_dot_estimation_EKF,'LineWidth',0.5, color = '[0.6350 0.0780 0.1840]');hold on;
plot(phi_dot_estimation_UKF,'LineWidth',0.5, color = 'g'); 
xlim([0,t_max/dt]); ylabel('phi\_dot');
legend({'phi\_dot reale','phi\_dot stimata EKF', 'phi\_dot stimata UKF'},'orientation','horizontal','location','southoutside');
title(t2,'EKF results'); xlabel(t2, 'Iterations');


figure(3);
t3 = tiledlayout(2,1);
nexttile; 
% Plot psi
plot(psi_real,'--','LineWidth',0.5, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(psi_estimation_EKF,'LineWidth',0.5, color = '[0.6350 0.0780 0.1840]');hold on;
plot(psi_estimation_UKF,'LineWidth',0.5, color = 'g'); 
xlim([0,t_max/dt]); ylabel('psi');
legend({'psi reale','psi stimata EKF','psi stimata UKF'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot psi_dot
plot(psi_dot_real,'--','LineWidth',0.5, color = '[0.3010 0.7450 0.9330]'); hold on; 
plot(psi_dot_estimation_EKF,'LineWidth',0.5, color = '[0.6350 0.0780 0.1840]');hold on;
plot(psi_dot_estimation_UKF,'LineWidth',0.5, color = 'g'); 
xlim([0,t_max/dt]); ylabel('psi_dot');
legend({'psi\_dot reale','psi\_dot stimata EKF','psi\_dot stimata UKF'},'orientation','horizontal','location','southoutside');
title(t3,'EKF results'); xlabel(t3, 'Iterations');

% Plot innovazione lungo tutte le componenti
figure(4);
t1 = tiledlayout(2,1);
nexttile; 
% Plot innovazione per psi
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(1,i) == 100
        plot(t,innovation_EKF(1,i),' ');hold on;
    else
        plot(t,innovation_EKF(1,i),'o','LineWidth',0.6, color = '[0.3010 0.7450 0.9330]');hold on;
    end
    if innovation_UKF(1,i) == 100
        plot(t,innovation_UKF(1,i),' ');hold on;
    else
        plot(t,innovation_UKF(1,i),'o','LineWidth',0.6, color = 'm');hold on;
    end
end
xlim([0,t_max]); 
legend('Innovazione per psi EKF','innovazione per psi UKF','orientation','horizontal','location','southoutside');

nexttile; 
% Plot innovazione per phi\_dot
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(2,i) == 100
        plot(t,innovation_EKF(2,i),' ');hold on;
    else
        plot(t,innovation_EKF(2,i),'o','LineWidth',0.6, color = '[0.3010 0.7450 0.9330]');hold on;
    end
    if innovation_UKF(2,i) == 100
        plot(t,innovation_UKF(2,i),' ');hold on;
    else
        plot(t,innovation_UKF(2,i),'o','LineWidth',0.6, color = 'm');hold on;
    end
end
xlim([0,t_max]); 
%plot(innovation_EKF(2,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on;
legend('Innovazione per phi\_dot EKF','Innovazione per phi\_dot UKF','orientation','horizontal','location','southoutside');
title(t1,'Innovazione'); xlabel(t1,'Time (s)');


figure(5);
t2 = tiledlayout(2,1);
nexttile; 
% Plot Innovazione per dx
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(3,i) == 100
        plot(t,innovation_EKF(3,i),' ');hold on;
    else
        plot(t,innovation_EKF(3,i),'o','LineWidth',0.6, color = '[0.3010 0.7450 0.9330]');hold on;
    end
    if innovation_UKF(3,i) == 100
        plot(t,innovation_UKF(3,i),' ');hold on;
    else
        plot(t,innovation_UKF(3,i),'o','LineWidth',0.6, color = 'm');hold on;
    end
end
xlim([0,t_max]); ylim([-1,1]);
%plot(innovation_EKF(3,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
legend({'Innovazione per dx EKF', 'Innovazione per dx UKF'},'orientation','horizontal','location','southoutside');

nexttile; 
% Plot innovazione per db
t = 0;
for i = 1 : size(innovation_EKF,2)
    t = t + dt;
    if innovation_EKF(4,i) == 100
        plot(t,innovation_EKF(4,i),' ');hold on;
    else
        plot(t,innovation_EKF(4,i),'o','LineWidth',0.6, color = '[0.3010 0.7450 0.9330]');hold on;
    end
    if innovation_UKF(4,i) == 100
        plot(t,innovation_UKF(4,i),' ');hold on;
    else
        plot(t,innovation_UKF(4,i),'o','LineWidth',0.6, color = 'm');hold on;
    end
end
xlim([0,t_max]); ylim([-1,1]);
%plot(innovation_EKF(4,:),'LineWidth',1, color = '[0.3010 0.7450 0.9330]'); hold on; 
legend('Innovazione per db EKF','Innovazione per db UKF','orientation','horizontal','location','southoutside');
title(t2,'Innovazione EKF'); xlabel(t2,'Time (s)');