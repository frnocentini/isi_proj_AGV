function F = function_F(in1)
%function_F
%    F = function_F(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    06-Dec-2022 18:16:03

phi_dot = in1(:,4);
psi_ = in1(:,5);
psi_dot = in1(:,6);
tau_phi_s = in1(:,7);
tau_psi_s = in1(:,8);
theta = in1(:,3);
w_tau_phi = in1(:,14);
w_tau_psi = in1(:,13);
t2 = cos(psi_);
t3 = cos(theta);
t4 = sin(psi_);
t5 = sin(theta);
t6 = psi_.*2.0;
t7 = cos(t6);
t8 = t2.^2;
t9 = t2.^3;
t10 = sin(t6);
t11 = t8.*5.357e+3;
t12 = t7.*2.6785e+4;
t13 = t8.*2.1428e+4;
t14 = t11+2.268e+3;
t15 = t13+9.072e+3;
t16 = 1.0./t14.^2;
t17 = 1.0./t15;
mt1 = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,phi_dot.*t2.*t5.*(-1.0./1.0e+2),(phi_dot.*t2.*t3)./1.0e+2,1.0,0.0,0.0,0.0,(t2.*t3)./1.0e+2,(t2.*t5)./1.0e+2,t4.*(-1.0./1.0e+2),(t12+psi_dot.*t10.*2.6785e+3+4.9465e+4)./(t12+4.9465e+4),0.0,psi_dot.*t2.*t17.*3.05e+2,phi_dot.*t3.*t4.*(-1.0./1.0e+2),phi_dot.*t4.*t5.*(-1.0./1.0e+2),phi_dot.*t2.*(-1.0./1.0e+2)];
mt2 = [(t16.*(phi_dot.*psi_dot.*-1.2149676e+7+t10.*tau_phi_s.*2.6785e+7+t2.*tau_psi_s.*6.491e+6-t9.*tau_psi_s.*2.6785e+6+t10.*w_tau_phi.*2.6785e+7+t2.*w_tau_psi.*6.491e+6-t9.*w_tau_psi.*2.6785e+6+phi_dot.*psi_dot.*t8.*5.2996801e+7))./1.0e+1,1.0,(t16.*(t2.*tau_phi_s.*2.5964e+6-t9.*tau_phi_s.*1.0714e+6+t10.*tau_psi_s.*1.525e+5+t2.*w_tau_phi.*2.5964e+6-t9.*w_tau_phi.*1.0714e+6+t10.*w_tau_psi.*1.525e+5-phi_dot.*psi_dot.*t4.^3.*1.633885e+6+phi_dot.*psi_dot.*t4.*9.42145e+5))./4.0,0.0,0.0,0.0];
mt3 = [(phi_dot.*t10.*5.357e+3)./(t7.*5.357e+4+9.893e+4),1.0./1.0e+1,t17.*(t15+phi_dot.*t2.*3.05e+2)];
F = reshape([mt1,mt2,mt3],6,6);