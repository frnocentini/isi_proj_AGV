function F = function_F(in1)
%function_F
%    F = function_F(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    08-Dec-2022 11:09:17

tau1 = in1(7,:);
tau2 = in1(8,:);
w1 = in1(13,:);
w2 = in1(14,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
t2 = cos(x3);
t3 = cos(x5);
t4 = sin(x3);
t5 = sin(x5);
t6 = x5.*2.0;
t7 = x5.*3.0;
t8 = cos(t6);
t9 = cos(t7);
t10 = t3.^2;
t11 = t3.^3;
t12 = sin(t6);
t13 = t10.*9.678517074195902e+18;
t14 = t10.*1.258207219645467e+21;
t15 = t8.*4.839258537097951e+18;
t16 = t13+2.009984661193124e+19;
t18 = t15+2.493910514902919e+19;
t19 = t14+2.612980059551061e+21;
t17 = 1.0./t16.^2;
t20 = 1.0./t18;
t21 = 1.0./t19;
et1 = t3.*tau2.*5.775184788405271e+36-t11.*tau2.*1.416615391387549e+36+t12.*tau1.*1.841600008803813e+37+t3.*w2.*5.775184788405271e+36;
et2 = t11.*w2.*-1.416615391387549e+36+t12.*w1.*1.841600008803813e+37-x4.*x6.*1.945367086222952e+37+t10.*x4.*x6.*4.827471100000921e+37;
et3 = t3.*tau1.*6.126540218323992e+38-t9.*tau1.*4.604000022009533e+37+t12.*tau2.*4.358569397017722e+37+t3.*w1.*6.126540218323992e+38;
et4 = t9.*w1.*-4.604000022009533e+37+t12.*w2.*4.358569397017722e+37+x4.*x6.*sin(t7).*7.205260034444927e+37+x4.*x6.*sin(x5.*5.0).*2.600496556373011e+21;
et5 = t5.*x4.*x6.*(-5.264879421010097e+38);
mt1 = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,t3.*t4.*x4.*(-1.0./1.0e+2),(t2.*t3.*x4)./1.0e+2,1.0,0.0,0.0,0.0,(t2.*t3)./1.0e+2,(t3.*t4)./1.0e+2,t5.*(-1.0./1.3e+2),t20.*(t18+t12.*x6.*4.839258537097953e+17),0.0,-t3.*t21.*x6.*(t10.*4.299e+3-2.977836368612715e+19),t2.*t5.*x4.*(-1.0./1.0e+2),t4.*t5.*x4.*(-1.0./1.0e+2),t3.*x4.*(-1.0./1.3e+2),t17.*(et1+et2),1.0,(t17.*(et3+et4+et5))./1.3e+2,0.0,0.0,0.0,t12.*t20.*x4.*4.839258537097953e+17,1.0./1.0e+1];
mt2 = [t21.*(t19+t3.*x4.*2.977836368612715e+19-t11.*x4.*4.299e+3+2.612980059551061e+21)];
F = reshape([mt1,mt2],6,6);
