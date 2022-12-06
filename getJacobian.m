function j_matrix = getJacobian()
    global rp l L IPy D IPz Mv mp ma IAy ra d IGz IAz a b
    syms x y theta phi_dot psi_ psi_dot %variabili di stato simboliche
    X = [x y theta phi_dot psi_ psi_dot]';%vettore di stato simbolico
    syms w_psi_ w_phi_dot w_dx w_db %rumore di misura simbolico
    %syms a rp b L l IPy IPz mp ma D IAz IGz d ra IAy %altre variabili simboliche utili nel modello
    syms tau_psi_s tau_phi_s %variabili di ingresso simboliche
    syms w_tau_psi w_tau_phi %rumore simbolico sulle variabili di ingresso
    
    load('dataset')
    dt = log_vars.dt;

    w_in = [w_tau_phi w_tau_psi]';
    w_sens = [w_psi_ w_phi_dot w_dx w_db]';

    %definizione della h (modello di osservazione)
    %vettore di variabili misurate [psi_,phi_dot,dx,db]
    h = [psi_+w_psi_,phi_dot+w_phi_dot,x+w_dx,sqrt((y^2)+(D-x)^2)+w_db]';
    
    %definizione della f (funzione di transizione di stato) e sua
    %discretizzazione
    %Matrice della dinamica
    B=[a*rp^2*cos(psi_)^2+b*(rp/L)^2*sin(psi_)^2+IPy -IPz*(rp/L)*sin(psi_);-IPz*(rp/L)*sin(psi_) IPz];
    C=[(rp^2*b/L^2-a*rp^2)*cos(psi_)*sin(psi_)*psi_dot*phi_dot;-IPz*rp/L*cos(psi_)*psi_dot*phi_dot];
    q_ddot = inv(B)*([tau_phi_s + w_tau_phi; tau_psi_s + w_tau_psi]-C);
    
    %X_dot
    x_dot = rp*phi_dot*cos(theta)*cos(psi_);
    y_dot = rp*phi_dot*sin(theta)*cos(psi_);
    theta_dot = -rp*phi_dot*sin(psi_)/L;
    phi_ddot = q_ddot(1);
    psi_dot = psi_dot;
    psi__ddot = q_ddot(2);
    f_cont = [x_dot y_dot theta_dot phi_ddot psi_dot psi__ddot]';
    f = X + dt*f_cont;   % funzione di stato f discretizzata
    F = simplify(jacobian(f', X'));
    T = simplify(jacobian(f', w_in'));
    H = simplify(jacobian(h', X'));
    M = simplify(jacobian(h', w_sens'));
    
    j_matrix = struct();
    j_matrix.F = F;
    j_matrix.T = T;
    j_matrix.H = H;
    j_matrix.M = M;
    j_matrix.f = f_cont;
    j_matrix.h = h;
end
