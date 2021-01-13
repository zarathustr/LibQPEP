% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Unified Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description 
%           in A Globally Optimal Framework.
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com
%
%
% test_cov_pTop.m: The globally optimal solution and covariance estimation
%                  of point-to-plane registration problem



clear all
close all
clc

if(verLessThan('matlab', '9.0.0'))
   error('The MATLAB version is too old to be supported.'); 
end

warning('off');
format long g

addpath('func_files');
addpath('solvers');
addpath('utils');
addpath('calib');
addpath('sdpt3');
run(fullfile('sdpt3', 'install_sdpt3.m'))
addpath('sedumi');
run(fullfile('sedumi', 'install_sedumi.m'))
p = genpath('YALMIP');
addpath(p);


m = 10;
[X, Y] = meshgrid(linspace(-2, 2, m), linspace(-2, 2, m));
X = reshape(X, 1, []);
Y = reshape(Y, 1, []);
Z = sin(X) .* cos(Y) .* cos(X);
len = m^2;

R0 = orthonormalize(randn(3, 3));
q0 = dcm2quat(R0).';
if(q0(1) < 0)
    q0 = - q0;
end
t0 = randn(3, 1);
X0 = [R0, t0;
     zeros(1, 3), 1];
 
r0 = [X; Y; Z].';
b0 = zeros(len, 3);
noise = 5e-100;
for i = 1 : len
    b0(i, :) = (R0 * r0(i, :).' + t0 + noise * randn(3, 1)).' + noise * randn(1, 3);
end

nvr = pcnormals(pointCloud(r0));
nvb = pcnormals(pointCloud(b0));

J = J_func_pTop(q0, t0, r0, b0, nvr)
J = J_func_pTop(q0, t0, r0, b0, nvb)


[W0, Q0, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = pTop_WQ_new(r0, b0, nvb);
t_funcs = {@t1_pTop_func_new, @t2_pTop_func_new, @t3_pTop_func_new};

W0_ = W0(1 : 3, :);
Q0_ = Q0(1 : 3, :);
W0_(1, :) = W0(1, :) + W0(2, :) + W0(3, :);
W0_(2, :) = W0(2, :) + W0(3, :) + W0(4, :);
W0_(3, :) = W0(3, :) + W0(4, :) + W0(1, :);
Q0_(1, :) = Q0(1, :) + Q0(2, :) + Q0(3, :);
Q0_(2, :) = Q0(2, :) + Q0(3, :) + Q0(4, :);
Q0_(3, :) = Q0(3, :) + Q0(4, :) + Q0(1, :);
X0

[~, ~, X] = QPEP_WQ_grobner(W0, Q0, @solver_WQ_approx, ...
                            @mon_J_pure_pTop_func_new, t_funcs, coef_J_pure, coefs_tq, pinvG, {[1, 2, 3]});
X_grobner = X

[~, ~, X_] = QPEP_lm_single(dcm2quat(X(1 : 3, 1 : 3)).', 100, 5e-2, ...
                          @eq_pTop_func_new, @Jacob_pTop_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
X_

num = 1000;
r_noise = 1e-3;
b_noise = 1e-3;
n_noise = 1e-3;
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

tic;
for j = 1 : num
    r = r0 + r_noise * randn(len, 3);
    b = b0 + b_noise * randn(len, 3);
    n = nvb + n_noise * randn(len, 3);

    [W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = pTop_WQDGc_new(r, b, n);
    Ws(:, :, j) = W;
    Qs(:, :, j) = Q;

    R = QPEP_WQ_grobner(W, Q, @solver_WQ_approx, ...
                        @mon_J_pure_pTop_func_new, t_funcs, coef_J_pure, coefs_tq, pinvG, {[1, 2, 3]});
    
    q = dcm2quat(real(R)).';
    if(q(1) < 0)
        q = - q;
    end
    [~, ~, X_] = QPEP_lm_single(q, 100, 5e-2, ...
                          @eq_pTop_func_new, @Jacob_pTop_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
    q = dcm2quat(X_(1 : 3, 1 : 3)).';
    t = X_(1 : 3, 4);
    if(q(1) < 0)
        q = - q;
    end
    
    qs(:, j) = q;
    ts(:, j) = t;
    y = y_func_pTop_new(q);
    v = v_func_pTop_new(q);
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
time = toc / num

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

r = r0 + r_noise * randn(len, 3);
b = b0 + b_noise * randn(len, 3);
n = nvb + n_noise * randn(len, 3);
    
[W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = pTop_WQDGc_new(r, b, n);
R = QPEP_WQ_grobner(W, Q, @solver_WQ_approx, ...
                    @mon_J_pure_pTop_func_new, t_funcs, coef_J_pure, coefs_tq, pinvG, {[1, 2, 3]});
R = real(R);
q = dcm2quat(R).';
if(q(1) < 0)
    q = - q;
end
    
[~, ~, X_] = QPEP_lm_single(q, 100, 5e-2, ...
                          @eq_pTop_func_new, @Jacob_pTop_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
q = dcm2quat(X_(1 : 3, 1 : 3)).';
if(q(1) < 0)
    q = - q;
end

y = y_func_pTop_new(q);
v = v_func_pTop_new(q);
V_cal = kron(v.', eye(3));
Y_cal = kron(y.', eye(3));

syms q0_ q1_ q2_ q3_
q_ = [q0_; q1_; q2_; q3_];
y_ = y_func_pTop_new(q_);
v_ = v_func_pTop_new(q_);

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
          
          
scalings = 1e3;
solver = 'sedumi';
verbose = true;
for i = 1 : length(scalings)
    scaling = scalings(i);

    epsX = scaling * 1e-23;
    epsQuat = scaling * 1e-15;
    
    
    [XX, ff] = optimize_quat_cov2(q, F, scaling * cov_left, solver, verbose);
    ff = ff / scaling^2
    scale = abs(norm(cov_left, 'fro') / norm(F * XX * F.', 'fro'));
    XX * scale
    Sigma_q_stat

    figure(i);
    plot_q_cov(Sigma_q_stat, XX * scale, qs);
    if(~ispc())
        set(gcf, 'Unit', 'Centimeters', 'Position', 6 * [5.5 5 5.5 5])
    end
end


          




