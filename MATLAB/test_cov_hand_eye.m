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
% test_cov_hand_eye.m: The globally optimal solution and covariance estimation
%                      of hand-eye calibration problem




clear all
close all
clc

warning('off');
format long g


addpath('func_files');
addpath('solvers');
addpath('utils');
addpath('calib');
addpath('sdpt3');
run(fullfile('sdpt3', 'install_sdpt3.m'))
addpath('sedumi');
if(~verLessThan('matlab', '8.0.0'))
    run(fullfile('sedumi', 'install_sedumi.m'))
end
p = genpath('YALMIP');
addpath(p);


len = 30;
A0 = zeros(4, 4, len);
B0 = zeros(4, 4, len);
R0 = orthonormalize(randn(3, 3));
q0 = positive_quat(dcm2quat(R0).');
t0 = randn(3, 1);
X0 = [R0, t0;
     zeros(1, 3), 1];
noise = 1e-3;
for i = 1 : len
    A0(1 : 3, 1 : 3, i) = orthonormalize(randn(3, 3));
    A0(1 : 3, 4, i) = randn(3, 1);
    A0(4, 4, i) = 1;
    B0(:, :, i) = inv(X0) * A0(:, :, i) * X0;
    B0(1 : 3, 1 : 3, i) = orthonormalize(B0(1 : 3, 1 : 3, i) + noise * randn(3, 3));
    B0(1 : 3, 4, i) = B0(1 : 3, 4, i) + noise * randn(3, 1);
end

[W0, Q0, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = hand_eye_WQ_new(A0, B0);
t_funcs = {@t1_hand_eye_func_new, @t2_hand_eye_func_new, @t3_hand_eye_func_new};

W0_ = W0(1 : 3, :);
Q0_ = Q0(1 : 3, :);
W0_(1, :) = W0(1, :) + W0(2, :) + W0(3, :);
W0_(2, :) = W0(2, :) + W0(3, :) + W0(4, :);
W0_(3, :) = W0(3, :) + W0(4, :) + W0(1, :);
Q0_(1, :) = Q0(1, :) + Q0(2, :) + Q0(3, :);
Q0_(2, :) = Q0(2, :) + Q0(3, :) + Q0(4, :);
Q0_(3, :) = Q0(3, :) + Q0(4, :) + Q0(1, :);
X0


[~, ~, X, ~, qs] = QPEP_WQ_grobner(W0, Q0, @solver_WQ, @mon_J_pure_hand_eye_func_new, ...
                            t_funcs, coef_J_pure, coefs_tq, pinvG, {[1, 2, 3]});
X_grobner = X
[~, ~, X_] = QPEP_lm_single(dcm2quat(X(1 : 3, 1 : 3)).', 1000, 5e-2, ...
                          @eq_hand_eye_func_new, @Jacob_hand_eye_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
X_
syms qq [4, 1]
uu = kron(qq, kron(qq, qq));
MM = jacobian(uu, qq);
u_func = matlabFunction(uu, 'Vars', {qq});
MM_func = matlabFunction(MM / 3, 'Vars', {qq});
iter = 1000;
q_trj = sym(zeros(4, iter));
% q = randn(4, 1); 
q = dcm2quat(X(1 : 3, 1 : 3)).';
q = q ./ norm(q);
vpa_len = 500;
q = vpa(q, vpa_len);
QQ0 = vpa(Q0, vpa_len);
for i = 1 : iter
    q_trj(:, iter) = q;
    M = vpa(W0 * MM_func(q), vpa_len);
    [V, D] = eig(M - QQ0);
    D
    q = V(:, 1);
    q = q ./ norm(q)
    QQ0 = QQ0 + D(1, 1) * eye(4);
end
q
dcm2quat(X(1 : 3, 1 : 3)).'

num = 5000;
quat_noise = 1e-2;
trans_noise = 1e-2;
ys = zeros(9, num);
vs = zeros(28, num);
Ds = zeros(3, 28, num);
ds = zeros(3 * 28, num);
Gs = zeros(3, 9, num);
gs = zeros(3 * 9, num);
cs = zeros(3, num);

Dvs = zeros(3, num);
Gys = zeros(3, num);
DvGys = zeros(3, num);

qs = zeros(4, num);
ts = zeros(3, num);
ps = zeros(3, num);
Ws = zeros(4, 64, num);
Qs = zeros(4, 4, num);
for j = 1 : num
    A = A0;
    B = B0;
    q = dcm2quat(A0(1 : 3, 1 : 3, :));
    term = [ones(len, 1), quat_noise * randn(len, 3)];
    term = norm_quat(term);
    q = quatmultiply(term, q);
    q = norm_quat(q);
    A(1 : 3, 1 : 3, :) = quat2dcm(q);
    
    q = dcm2quat(B0(1 : 3, 1 : 3, :));
    term = [ones(len, 1), quat_noise * randn(len, 3)];
    term = norm_quat(term);
    q = quatmultiply(term, q);
    q = norm_quat(q);
    B(1 : 3, 1 : 3, :) = quat2dcm(q);
    
    A(1 : 3, 4, :) = A0(1 : 3, 4, :) + trans_noise * randn(3, 1, len);
    B(1 : 3, 4, :) = B0(1 : 3, 4, :) + trans_noise * randn(3, 1, len);
    
    [W, Q] = hand_eye_WQ_new(A, B);
    Ws(:, :, j) = W;
    Qs(:, :, j) = Q;
    [D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = hand_eye_DGc_new(A, B);
    R = QPEP_WQ_grobner(W, Q, @solver_WQ_approx, @mon_J_pure_hand_eye_func_new, ...
                        t_funcs, coef_J_pure, coefs_tq, pinvG, {[1, 2, 3]});
    
    q = positive_quat(dcm2quat(real(R)).');
    [~, ~, X_] = QPEP_lm_single(q, 1000, 5e-2, ...
                          @eq_hand_eye_func_new, @Jacob_hand_eye_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
    q = positive_quat(dcm2quat(X_(1 : 3, 1 : 3)).');
    t = X_(1 : 3, 4);
    
    qs(:, j) = q;
    ts(:, j) = t;
    y = y_func_hand_eye_new(q);
    v = v_func_hand_eye_new(q);
    ys(:, j) = y;
    vs(:, j) = v;
    Ds(:, :, j) = D;
    ds(:, j) = vec(D);
    Gs(:, :, j) = G;
    gs(:, j) = vec(G);
    cs(:, j) = c;
    Dvs(:, j) = D * v;
    Gys(:, j) = G * y;
    ps(:, j) = D * v + G * y + c;
    DvGys(:, j) = D * v + G * y;
    
    if(mod(j, 1000) == 0)
        fprintf('%d\n', j);
    end
end

Sigma_q_stat = covx(qs.', qs.');
Sigma_y_stat = covx(ys.', ys.');
Sigma_v_stat = covx(vs.', vs.');

Sigma_d_stat = covx(ds.', ds.');
Sigma_g_stat = covx(gs.', gs.');
Sigma_c_stat = covx(cs.', cs.');

Sigma_d_q_stat = covx(ds.', qs.');
Sigma_q_d_stat = covx(qs.', ds.');
Sigma_g_q_stat = covx(gs.', qs.');
Sigma_q_g_stat = covx(qs.', gs.');
Sigma_c_q_stat = covx(cs.', qs.');
Sigma_q_c_stat = covx(qs.', cs.');

Sigma_d_g_stat = covx(ds.', gs.');
Sigma_g_d_stat = covx(gs.', ds.');
Sigma_d_c_stat = covx(ds.', cs.');
Sigma_c_d_stat = covx(cs.', ds.');
Sigma_c_g_stat = covx(cs.', gs.');
Sigma_g_c_stat = covx(gs.', cs.');

Sigma_p_stat = covx(ps.', ps.');

Sigma_Dv_stat = covx(Dvs.', Dvs.');
Sigma_Gy_stat = covx(Gys.', Gys.');
Sigma_DvGy_stat = covx(DvGys.', DvGys.');

mean_D = zeros(3, 28);
mean_d = mean(ds, 2);
mean_G = zeros(3, 9);
mean_g = mean(gs, 2);
mean_Dv = mean(Dvs, 2);
mean_Gy = mean(Gys, 2);
mean_W = zeros(4, 64);
mean_Q = zeros(4, 4);

for i = 1 : num
    mean_D = mean_D + 1 / num * Ds(:, :, i);
    mean_G = mean_G + 1 / num * Gs(:, :, i);
    mean_W = mean_W + 1 / num * Ws(:, :, i);
    mean_Q = mean_Q + 1 / num * Qs(:, :, i);
end

len_ = len;

A = A0;
B = B0;
q = dcm2quat(A0(1 : 3, 1 : 3, :));
term = [ones(len, 1), quat_noise * randn(len, 3)];
term = norm_quat(term);
q = quatmultiply(term, q);
q = norm_quat(q);
A(1 : 3, 1 : 3, :) = quat2dcm(q);
    
q = dcm2quat(B0(1 : 3, 1 : 3, :));
term = [ones(len, 1), quat_noise * randn(len, 3)];
term = norm_quat(term);
q = quatmultiply(term, q);
q = norm_quat(q);
B(1 : 3, 1 : 3, :) = quat2dcm(q);
    
[W, Q] = hand_eye_WQ_new(A, B);
[D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = hand_eye_DGc_new(A, B);
R = QPEP_WQ_grobner(W, Q, @solver_WQ_approx, @mon_J_pure_hand_eye_func_new, ...
                    t_funcs, coef_J_pure, coefs_tq, pinvG, {[1, 2, 3]});
R = real(R);
q = positive_quat(dcm2quat(R).');
    
[~, ~, X_] = QPEP_lm_single(q, 1000, 5e-2, ...
                          @eq_hand_eye_func_new, @Jacob_hand_eye_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
q = positive_quat(dcm2quat(X_(1 : 3, 1 : 3)).');

y = y_func_hand_eye_new(q);
v = v_func_hand_eye_new(q);
V_cal = kron(v.', eye(3));
Y_cal = kron(y.', eye(3));

syms q0_ q1_ q2_ q3_
q_ = [q0_; q1_; q2_; q3_];
y_ = y_func_hand_eye_new(q_);
v_ = v_func_hand_eye_new(q_);

partial_v_q_sym = jacobian(v_, q_);
partial_v_q_sym_func = matlabFunction(partial_v_q_sym, 'Vars', {q_});
partial_v_q_sym_val = partial_v_q_sym_func(q);

partial_y_q_sym = jacobian(y_, q_);
partial_y_q_sym_func = matlabFunction(partial_y_q_sym, 'Vars', {q_});
partial_y_q_sym_val = partial_y_q_sym_func(q);

F = D * partial_v_q_sym_val + G * partial_y_q_sym_val;

cov_left = V_cal * Sigma_d_stat * V_cal.' + Y_cal * Sigma_g_stat * Y_cal.' + Sigma_c_stat + ...
           V_cal * Sigma_d_g_stat * Y_cal.' + V_cal * Sigma_d_c_stat + Y_cal * Sigma_g_c_stat + ...
           Sigma_c_d_stat * V_cal.' + Sigma_c_g_stat * Y_cal.' + Y_cal * Sigma_g_d_stat * V_cal.';

scalings = 1;
solver = 'sedumi';
verbose = true;
for i = 1 : length(scalings)
    scaling = scalings(i);
    
    [XX, ff] = optimize_quat_cov2(q, F, scaling * cov_left, solver, verbose);
    ff = ff / scaling^2
    scale = abs(norm(cov_left, 'fro') / norm(F * XX * F.', 'fro'));
    XX * scale
    Sigma_q_stat

    figure(i);
    plot_q_cov(Sigma_q_stat, XX * scale, qs);
end


          
