function [W, Q, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = hand_eye_WQ_new(A, B)
len = size(A, 3);
coef_J_pure = zeros(1, 85);
for i = 1 : len
    AA = A(1 : 3, :, i);
    BB = B(1 : 3, :, i);
    coef_J_pure = coef_J_pure + 1 / len * coef_J_pure_hand_eye_func_new(AA, BB);
end

G = G_hand_eye_func_new(coef_J_pure);
coeftq1 = coeftq1_hand_eye_func_new(coef_J_pure);
coeftq2 = coeftq2_hand_eye_func_new(coef_J_pure);
coeftq3 = coeftq3_hand_eye_func_new(coef_J_pure);
coef_Jacob1_qt = coef_Jacob1_qt_hand_eye_func_new(coef_J_pure);
coef_Jacob2_qt = coef_Jacob2_qt_hand_eye_func_new(coef_J_pure);
coef_Jacob3_qt = coef_Jacob3_qt_hand_eye_func_new(coef_J_pure);
coef_Jacob4_qt = coef_Jacob4_qt_hand_eye_func_new(coef_J_pure);


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

W = W_hand_eye_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
Q = Q_hand_eye_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);
coef_f_q_sym = coef_f_q_sym_hand_eye_func_new(pinvG, coefs_tq, coef_Jacob_qt_syms);          
end