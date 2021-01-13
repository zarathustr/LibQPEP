function [W, Q, D, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = pnp_WQD_new(image_pt, world_pt, K, scale)
len = size(image_pt, 1);
fx = K(1, 1); fy = K(2, 2);
cx = K(3, 1); cy = K(3, 2);
coef_J_pure = zeros(1, 70);
for i = 1 : len
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates
    image_point = image_pt(i, :);
   
    coef_J_pure = coef_J_pure + 1 / len * coef_J_pure_pnp_func_new(image_point, world_point, [fx, fy, cx, cy], scale);
end

G = G_pnp_func_new(coef_J_pure);
coeftq1 = coeftq1_pnp_func_new(coef_J_pure);
coeftq2 = coeftq2_pnp_func_new(coef_J_pure);
coeftq3 = coeftq3_pnp_func_new(coef_J_pure);
coef_Jacob1_qt = coef_Jacob1_qt_pnp_func_new(coef_J_pure);
coef_Jacob2_qt = coef_Jacob2_qt_pnp_func_new(coef_J_pure);
coef_Jacob3_qt = coef_Jacob3_qt_pnp_func_new(coef_J_pure);
coef_Jacob4_qt = coef_Jacob4_qt_pnp_func_new(coef_J_pure);

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

W = W_pnp_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
Q = Q_pnp_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
coef_f_q_sym = coef_f_q_sym_pnp_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
D = D_func_pnp_new(coef_f_q_sym(1, :), coef_f_q_sym(2, :), coef_f_q_sym(3, :), coef_f_q_sym(4, :));
end