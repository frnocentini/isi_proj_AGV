function T = function_T(in1)
%function_T
%    T = function_T(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    06-Dec-2022 17:54:32

psi_ = in1(:,5);
t2 = sin(psi_);
t3 = t2.^2;
t4 = t3.*5.357e+3;
t5 = t4-7.625e+3;
t6 = 1.0./t5;
t7 = t2.*t6.*5.0e+1;
t8 = -t7;
T = reshape([0.0,0.0,0.0,t6.*-5.0e+2,0.0,t8,0.0,0.0,0.0,t8,0.0,t6.*(t3.*5.356e+3-7.625e+3).*5.0],[6,2]);
