function H = function_H(in1)
%function_H
%    H = function_H(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    06-Dec-2022 18:16:04

x = in1(:,1);
y = in1(:,2);
t2 = y.^2;
t3 = x-2.1e+1./5.0;
t4 = t3.^2;
t5 = t2+t4;
t6 = 1.0./sqrt(t5);
H = reshape([0.0,0.0,1.0,(t6.*(x.*5.0-2.1e+1))./5.0,0.0,0.0,0.0,t6.*y,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,6]);
