function h = observation_model_h(in1)
%OBSERVATION_MODEL_H
%    H = OBSERVATION_MODEL_H(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    06-Dec-2022 18:16:03

phi_dot = in1(:,4);
psi_ = in1(:,5);
w_db = in1(:,10);
w_dx = in1(:,9);
w_psi_ = in1(:,7);
w_phi_dot = in1(:,8);
x = in1(:,1);
y = in1(:,2);
h = [conj(psi_)+conj(w_psi_);conj(phi_dot)+conj(w_phi_dot);conj(w_dx)+conj(x);conj(sqrt((x-2.1e+1./5.0).^2+y.^2))+conj(w_db)];
