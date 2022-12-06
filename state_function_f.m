function f = state_function_f(in1)
%STATE_FUNCTION_F
%    F = STATE_FUNCTION_F(IN1)

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
t2 = conj(phi_dot);
t3 = conj(psi_);
t4 = conj(psi_dot);
t5 = conj(tau_phi_s);
t6 = conj(tau_psi_s);
t7 = conj(theta);
t8 = conj(w_tau_phi);
t9 = conj(w_tau_psi);
t10 = cos(t3);
t11 = sin(t3);
t12 = t10.^2;
t13 = t11.^2;
t16 = (t2.*t4.*t10)./5.0e+2;
t18 = t2.*t4.*t10.*t11.*1.0712;
t14 = t13.*2.168e+3;
t15 = t12.*7.525e+3;
t17 = t6+t9+t16;
t20 = t5+t8+t18;
t19 = t14+t15+1.0e+2;
t21 = 1.0./t19;
f = [(t2.*t10.*cos(t7))./1.0e+1;(t2.*t10.*sin(t7))./1.0e+1;t2.*t11.*(-1.0./1.0e+1);t20.*t21.*5.0e+3+t11.*t17.*t21.*5.0e+2;t4;t11.*t20.*t21.*5.0e+2+t17.*t21.*(t13.*2.169e+3+t15+1.0e+2).*5.0e+1];