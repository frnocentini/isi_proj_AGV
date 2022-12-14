function F = jacobian_F(in1)
%jacobian_F
%    F = jacobian_F(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    14-Dec-2022 11:35:21

tau1 = in1(7,:);
tau2 = in1(8,:);
w1 = in1(13,:);
w2 = in1(14,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
et1 = tau1.*cos(x5.*3.0).*-4.604000022009533e+37-w1.*cos(x5.*3.0).*4.604000022009533e+37+tau2.*sin(x5.*2.0).*4.358569397017722e+37+w2.*sin(x5.*2.0).*4.358569397017722e+37;
et2 = tau1.*cos(x5).*6.126540218323992e+38+w1.*cos(x5).*6.126540218323992e+38-x4.*x6.*sin(x5).*5.264879421010097e+38+x4.*x6.*sin(x5.*3.0).*7.205260034444927e+37;
et3 = x4.*x6.*sin(x5.*5.0).*2.600496556373011e+21;
et4 = tau2.*cos(x5).^3.*-1.416615391387549e+36-w2.*cos(x5).^3.*1.416615391387549e+36+tau1.*sin(x5.*2.0).*1.841600008803813e+37+w1.*sin(x5.*2.0).*1.841600008803813e+37;
et5 = x4.*x6.*-1.945367086222952e+37+tau2.*cos(x5).*5.775184788405271e+36+w2.*cos(x5).*5.775184788405271e+36+x4.*x6.*cos(x5).^2.*4.827471100000921e+37;
mt1 = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,x4.*cos(x5).*sin(x3).*(-1.0./1.0e+2),(x4.*cos(x3).*cos(x5))./1.0e+2,1.0,0.0,0.0,0.0,(cos(x3).*cos(x5))./1.0e+2,(cos(x5).*sin(x3))./1.0e+2,sin(x5).*(-1.0./1.3e+2),(cos(x5.*2.0).*4.839258537097951e+18+x6.*sin(x5.*2.0).*4.839258537097953e+17+2.493910514902919e+19)./(cos(x5.*2.0).*4.839258537097951e+18+2.493910514902919e+19),0.0];
mt2 = [-(x6.*cos(x5).*(cos(x5).^2.*4.299e+3-2.977836368612715e+19))./(cos(x5).^2.*1.258207219645467e+21+2.612980059551061e+21),x4.*cos(x3).*sin(x5).*(-1.0./1.0e+2),x4.*sin(x3).*sin(x5).*(-1.0./1.0e+2),x4.*cos(x5).*(-1.0./1.3e+2),(et4+et5).*1.0./(cos(x5).^2.*9.678517074195902e+18+2.009984661193124e+19).^2,1.0];
mt3 = [(1.0./(cos(x5).^2.*9.678517074195902e+18+2.009984661193124e+19).^2.*(et1+et2+et3))./1.3e+2,0.0,0.0,0.0,(x4.*sin(x5.*2.0).*4.839258537097953e+17)./(cos(x5.*2.0).*4.839258537097951e+18+2.493910514902919e+19),1.0./1.0e+1];
mt4 = [(x4.*cos(x5).^3.*-4.299e+3+cos(x5).^2.*1.258207219645467e+21+x4.*cos(x5).*2.977836368612715e+19+2.612980059551061e+21)./(cos(x5).^2.*1.258207219645467e+21+2.612980059551061e+21)];
F = reshape([mt1,mt2,mt3,mt4],6,6);
