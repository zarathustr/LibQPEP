% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Globally Optimal Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description.
%           IEEE Transactions on Robotics.
%           https://doi.org/10.1109/TRO.2022.3155880
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com
%
%
% syms_pnp.m: Generating symbolic expressions for solving PnP problem




clear all
close all
clc

addpath('func_files');
addpath('utils');

syms b1 b2 real
syms r1 r2 r3 real

image_pt = [b1, b2];
world_pt = [r1, r2, r3];

syms fx fy cx cy real
K = [fx, 0, 0;
      0, fy, 0;
     cx, cy, 1];
  
syms q0 q1 q2 q3 real
syms t1 t2 t3 lambda scale real
q = [q0; q1; q2; q3];
tt = [t1; t2; t3];
RR = q2R(q);
XX = [
    RR, tt;
    zeros(1, 3), 1];
J_pure = J_pnp_loss(image_pt, world_pt, K, RR, tt) * scale;
[coef_J_pure, mon_J_pure] = coeffs(J_pure, [q; tt]);
generateFuncFile(coef_J_pure, fullfile('func_files', 'coef_J_pure_pnp_func_new.m'), {image_pt, world_pt, [fx, fy, cx, cy], scale});
generateFuncFile(mon_J_pure, fullfile('func_files', 'mon_J_pure_pnp_func_new.m'), {q, tt});
J = J_pure - 0.5 * lambda * (q.' * q - 1);
x = [q; tt; lambda];
Jacob = jacobian(J, x);
Jacob = expand(Jacob.');

coef_len = length(coef_J_pure);
str = sprintf('coef_J = sym(''coef_J'', [1, %d]);', coef_len);
eval(str);
showSyms(symvar(coef_J));
J = coef_J * mon_J_pure.' - 0.5 * lambda * (q.' * q - 1);
x = [q; tt; lambda];
Jacob = jacobian(J, x);
Jacob = expand(Jacob.');

ts = solveSyms(Jacob(5 : 7), tt);
ss = Jacob(1 : 4);
ss = subs(ss, q0^3, (q0 * (1 - q1^2 - q2^2 - q3^2)));
ss = subs(ss, q0^2, (1 - q1^2 - q2^2 - q3^2));
Jacob(1 : 4) = ss;

ss = subs(ss, t1, ts.t1);
ss = subs(ss, t2, ts.t2);
ss = subs(ss, t3, ts.t3);


[coeft1, mont1] = coeffs(Jacob(5), tt);
g1_ = coeft1(1 : 3).';
g1 = vpa(zeros(3, 1));
for i = 1 : 3
    idx = sscanf(char(mont1(i)), 't%d');
    if(~isempty(idx))
        g1(idx) = g1_(i);
    end
end
[coeft2, mont2] = coeffs(Jacob(6), tt);
g2_ = coeft2(1 : 3).';
g2 = vpa(zeros(3, 1));
for i = 1 : 3
    idx = sscanf(char(mont2(i)), 't%d');
    if(~isempty(idx))
        g2(idx) = g2_(i);
    end
end
[coeft3, mont3] = coeffs(Jacob(7), tt);
g3_ = coeft3(1 : 3).';
g3 = vpa(zeros(3, 1));
for i = 1 : 3
    idx = sscanf(char(mont3(i)), 't%d');
    if(~isempty(idx))
        g3(idx) = g3_(i);
    end
end

[coeftq1, montq1] = coeffs(expand(Jacob(5) - g1.' * tt ), q);
generateFuncFile(coeftq1, fullfile('func_files', 'coeftq1_pnp_func_new.m'), {coef_J});
[coeftq2, montq2] = coeffs(expand(Jacob(6) - g2.' * tt ), q);
generateFuncFile(coeftq2, fullfile('func_files', 'coeftq2_pnp_func_new.m'), {coef_J});
[coeftq3, montq3] = coeffs(expand(Jacob(7) - g3.' * tt ), q);
generateFuncFile(coeftq3, fullfile('func_files', 'coeftq3_pnp_func_new.m'), {coef_J});
G = - [g1.'; g2.'; g3.'];
generateFuncFile(G, fullfile('func_files', 'G_pnp_func_new.m'), {coef_J});

pinvG = sym('pinvG', [3, 3]);
showSyms(symvar(pinvG));
coefs_tq = sym('coefs_tq', [3, 10]);
showSyms(symvar(coefs_tq));
ts = pinvG * coefs_tq * montq1.';
t1 = ts(1);
t2 = ts(2);
t3 = ts(3);
generateFuncFile(t1, fullfile('func_files', 't1_pnp_func_new.m'), {pinvG, coefs_tq, q});
generateFuncFile(t2, fullfile('func_files', 't2_pnp_func_new.m'), {pinvG, coefs_tq, q});
generateFuncFile(t3, fullfile('func_files', 't3_pnp_func_new.m'), {pinvG, coefs_tq, q});

[coef_Jacob1_qt, mon_Jacob1_qt] = coeffs(Jacob(1) + lambda * q0, [q; tt]);
generateFuncFile(coef_Jacob1_qt, fullfile('func_files', 'coef_Jacob1_qt_pnp_func_new.m'), {coef_J});
[coef_Jacob2_qt, mon_Jacob2_qt] = coeffs(Jacob(2) + lambda * q1, [q; tt]);
generateFuncFile(coef_Jacob2_qt, fullfile('func_files', 'coef_Jacob2_qt_pnp_func_new.m'), {coef_J});
[coef_Jacob3_qt, mon_Jacob3_qt] = coeffs(Jacob(3) + lambda * q2, [q; tt]);
generateFuncFile(coef_Jacob3_qt, fullfile('func_files', 'coef_Jacob3_qt_pnp_func_new.m'), {coef_J});
[coef_Jacob4_qt, mon_Jacob4_qt] = coeffs(Jacob(4) + lambda * q3, [q; tt]);
generateFuncFile(coef_Jacob4_qt, fullfile('func_files', 'coef_Jacob4_qt_pnp_func_new.m'), {coef_J});

coef_Jacob1_qt_syms = sym('coef_Jacob1_qt_syms', [1, 32]);
coef_Jacob2_qt_syms = sym('coef_Jacob2_qt_syms', [1, 32]);
coef_Jacob3_qt_syms = sym('coef_Jacob3_qt_syms', [1, 32]);
coef_Jacob4_qt_syms = sym('coef_Jacob4_qt_syms', [1, 32]);
showSyms(symvar(coef_Jacob1_qt_syms));
showSyms(symvar(coef_Jacob2_qt_syms));
showSyms(symvar(coef_Jacob3_qt_syms));
showSyms(symvar(coef_Jacob4_qt_syms));
coef_Jacob_qt_syms = [
    coef_Jacob1_qt_syms;
    coef_Jacob2_qt_syms;
    coef_Jacob3_qt_syms;
    coef_Jacob4_qt_syms;
    ];
Jacob_ = coef_Jacob_qt_syms * mon_Jacob1_qt.';
    
eqs = expand(eval([Jacob_; Jacob(8)]));

eqs_ = eqs;
eqs_(1 : 4) = eval(subs(eqs_(1 : 4), q0^2, (1 - q3 * q3 - q1 * q1 - q2 * q2)));
eqs_(1 : 4) = eval(subs(eqs_(1 : 4), q0^3, q0 * (1 - q3 * q3 - q1 * q1 - q2 * q2)));

f0 = eqs(1);
f1 = eqs(2);
f2 = eqs(3);
f3 = eqs(4);
[coef_f0_q, mon_f0_q] = coeffs(f0, q);
[coef_f1_q, mon_f1_q] = coeffs(f1, q);
[coef_f2_q, mon_f2_q] = coeffs(f2, q);
[coef_f3_q, mon_f3_q] = coeffs(f3, q);
coef_f0_q_sym = sym('coef_f0_q_sym', [1, 24]);
coef_f1_q_sym = sym('coef_f1_q_sym', [1, 24]);
coef_f2_q_sym = sym('coef_f2_q_sym', [1, 24]);
coef_f3_q_sym = sym('coef_f3_q_sym', [1, 24]);
showSyms(symvar(coef_f0_q_sym));
showSyms(symvar(coef_f1_q_sym));
showSyms(symvar(coef_f2_q_sym));
showSyms(symvar(coef_f3_q_sym));
f0 = coef_f0_q_sym * mon_f0_q.';
f1 = coef_f1_q_sym * mon_f1_q.';
f2 = coef_f2_q_sym * mon_f2_q.';
f3 = coef_f3_q_sym * mon_f3_q.';
coef_f_q_sym = [
    coef_f0_q;
    coef_f1_q;
    coef_f2_q;
    coef_f3_q;
    ];
generateFuncFile(coef_f_q_sym, fullfile('func_files', 'coef_f_q_sym_pnp_func_new.m'), {pinvG, coefs_tq, coef_Jacob_qt_syms});

eq = expand([
     f0 * q1 - f1 * q0;
     f0 * q2 - f2 * q0;
     f0 * q3 - f3 * q0;
     q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3 - 1;
     ]);
generateFuncFile(eq, fullfile('func_files', 'eq_pnp_func_new.m'), {coef_f0_q_sym, coef_f1_q_sym, coef_f2_q_sym, coef_f3_q_sym, q});


eq_ = eq;
[coefeq1, moneq1] = coeffs(eq_(1), q); 
[coefeq2, moneq2] = coeffs(eq_(2), q); 
[coefeq3, moneq3] = coeffs(eq_(3), q); 

v1 = moneq1.';
generateFuncFile(v1, fullfile('func_files', 'v1_func_pnp_new.m'), {q});
v2 = moneq2.';
generateFuncFile(v2, fullfile('func_files', 'v2_func_pnp_new.m'), {q});
v3 = moneq3.';
generateFuncFile(v3, fullfile('func_files', 'v3_func_pnp_new.m'), {q});

D = [
    coefeq1;
    coefeq2;
    coefeq3;
    ];
generateFuncFile(D, fullfile('func_files', 'D_func_pnp_new.m'), {coef_f0_q_sym, coef_f1_q_sym, coef_f2_q_sym, coef_f3_q_sym});

HH = jacobian(eq, q);
generateFuncFile(HH, fullfile('func_files', 'Jacob_pnp_func_new.m'), {coef_f0_q_sym, coef_f1_q_sym, coef_f2_q_sym, coef_f3_q_sym, q});
gradient = expand(HH.' * eq);

[coef1, mon1] = coeffs(eqs_(1), x);
[coef2, mon2] = coeffs(eqs_(2), x);
[coef3, mon3] = coeffs(eqs_(3), x);
[coef4, mon4] = coeffs(eqs_(4), x);

Q_sym = [
    coef1(7), coef1(14), coef1(18), coef1(20);
    coef2(7), coef2(14), coef2(18), coef2(20);
    coef3(7), coef3(14), coef3(18), coef3(20);
    coef4(7), coef4(14), coef4(18), coef4(20);
    ];
generateFuncFile(Q_sym, fullfile('func_files', 'Q_pnp_func_new.m'), {pinvG, coefs_tq, coef_Jacob_qt_syms});

res = expand(eqs_(1 : 4) - Q_sym * q);
H = jacobian(res, q);
h1 = H(:, 1);
h2 = H(:, 2);
h3 = H(:, 3);
h4 = H(:, 4);
P1 = jacobian(h1, q);
P2 = jacobian(h2, q);
P3 = jacobian(h3, q);
P4 = jacobian(h4, q);

for i = 1 : 4
    for j = 1 : 4
        str = sprintf('W%d%d = [jacobian(P%d(1, :), q%d).''; jacobian(P%d(2, :), q%d).''; jacobian(P%d(3, :), q%d).''; jacobian(P%d(4, :), q%d).''];', ...
                        i, j, i, j - 1, i, j - 1, i, j - 1, i, j - 1);
        eval(str);
    end
end


W1 = [W11, W12, W13, W14];
W2 = [W21, W22, W23, W24];
W3 = [W31, W32, W33, W34];
W4 = [W41, W42, W43, W44];
v = kron(q, q);
u = kron(v, q);
W = - [W1, W2, W3, W4] / 6;
generateFuncFile(W, fullfile('func_files', 'W_pnp_func_new.m'), {pinvG, coefs_tq, coef_Jacob_qt_syms});





