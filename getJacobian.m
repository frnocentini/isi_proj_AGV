function j_matrix = getJacobian()
    global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b
%     syms x y theta X(4) X(5) X(5)dot %variabili di stato simboliche
%     X = [x y theta X(4) X(5) X(5)dot]';%vettore di stato simbolico
%     syms V(1) V(2) w_dx w_db %rumore di misura simbolico
%     %syms a rp b L l IPy IPz mp ma D IAz IGz d ra IAy %altre variabili simboliche utili nel modello
%     syms tau_X(5)s tau_phi_s %variabili di ingresso simboliche
%     syms w_tau_psi w_tau_phi %rumore simbolico sulle variabili di ingresso

%% Vettore di stato X = [x y theta phi_dot psi psi_dot]
%Tau=[tau_phi tau_psi]  ingressi
Tau=sym('tau',[2 1],'real');

%W = [w_tau_phi w_tau_psi]  disturbi di processo
W = sym('w',[2 1],'real');

%V=[v_psi v_phi_dot v_dx v_db] disturbi di misura
V = sym('v',[4 1],'real');

%X vettore di stato
X = sym('x',[6 1],'real');

load('dataset')
dt = log_vars.dt;

%     w_in = [w_tau_phi w_tau_psi]';
%     w_sens = [V(1) V(2) w_dx w_db]';

%definizione della h (modello di osservazione)
%vettore di variabili misurate [X(5),X(4),dx,db]
h = [X(5)+V(1),X(4)+V(2),X(1)+V(3),sqrt((X(2)^2)+(D-X(1))^2)+V(4)]';

%definizione della f (funzione di transizione di stato) e sua
%discretizzazione
%Matrice della dinamica
B=[a*rp^2*cos(X(5))^2+b*(rp/L)^2*sin(X(5))^2+IPy -IPz*(rp/L)*sin(X(5));-IPz*(rp/L)*sin(X(5)) IPz];
C=[(rp^2*b/L^2-a*rp^2)*cos(X(5))*sin(X(5))*X(6)*X(4);-IPz*rp/L*cos(X(5))*X(6)*X(4)];
q_ddot = inv(B)*([Tau(1) + W(1); Tau(2) + W(2)]-C);
    
%X_dot
x1_dot = rp*X(4)*cos(X(3))*cos(X(5));
x2_dot = rp*X(4)*sin(X(3))*cos(X(5));
x3_dot = -rp*X(4)*sin(X(5))/L;
x4_dot = q_ddot(1);
x5_dot = X(6);
x6_dot = q_ddot(2);
f_cont = [x1_dot x2_dot x3_dot x4_dot x5_dot x6_dot]';
f_disc = X + dt*f_cont;
F = simplify(jacobian(f_disc, X));
T = simplify(jacobian(f_disc, W));
H = simplify(jacobian(h, X));
M = simplify(jacobian(h, V));

j_matrix = struct();
j_matrix.F = F;
j_matrix.T = T;
j_matrix.H = H;
j_matrix.M = M;
j_matrix.f = f_cont;
j_matrix.h = h;

end
