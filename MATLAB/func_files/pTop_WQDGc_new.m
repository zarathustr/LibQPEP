function [W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = pTop_WQDGc_new(r, b, nv)
len = size(r, 1);
coef_J_pure = zeros(1, 85);
for i = 1 : len
    rr = r(i, :).';
    bb = b(i, :).';
    nn = nv(i, :).';
   
    coef_J_pure = coef_J_pure + 1 / len * coef_J_pure_pTop_func_new(rr, bb, nn);
end
G = G_pTop_func_new(coef_J_pure);
coeftq1 = coeftq1_pTop_func_new(coef_J_pure);
coeftq2 = coeftq2_pTop_func_new(coef_J_pure);
coeftq3 = coeftq3_pTop_func_new(coef_J_pure);
coef_Jacob1_qt = coef_Jacob1_qt_pTop_func_new(coef_J_pure);
coef_Jacob2_qt = coef_Jacob2_qt_pTop_func_new(coef_J_pure);
coef_Jacob3_qt = coef_Jacob3_qt_pTop_func_new(coef_J_pure);
coef_Jacob4_qt = coef_Jacob4_qt_pTop_func_new(coef_J_pure);

pinvG = pinv(G); 
coefs_tq = [
    coeftq1;
    coeftq2;
    coeftq3;
    ];
coef_Jacob_qt_syms = [
    coef_Jacob1_qt;
    coef_Jacob2_qt;
    coef_Jacob3_qt;
    coef_Jacob4_qt;
    ];

coef_f_q_sym = coef_f_q_sym_pTop_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
W = W_pTop_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
Q = Q_pTop_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
D = D_func_pTop_new(coef_f_q_sym(1, :), coef_f_q_sym(2, :), coef_f_q_sym(3, :), coef_f_q_sym(4, :));
G = G_func_pTop_new(coef_f_q_sym(1, :), coef_f_q_sym(2, :), coef_f_q_sym(3, :), coef_f_q_sym(4, :));
c = c_func_pTop_new(coef_f_q_sym(1, :), coef_f_q_sym(2, :), coef_f_q_sym(3, :), coef_f_q_sym(4, :));
end