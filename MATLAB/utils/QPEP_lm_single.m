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



function [R, t, X] = QPEP_lm_single(q0, max_iter, mu, ...
                                  eq_func, Jacob_func, t_funcs, ...
                                  coef_f_q_sym, coefs_tq, pinvG)
eq_fun = @(q)(eq_func(coef_f_q_sym(1, :), ...
                               coef_f_q_sym(2, :), ...
                               coef_f_q_sym(3, :), ...
                               coef_f_q_sym(4, :), ...
                               q));
Jacob_fun = @(q)(Jacob_func(coef_f_q_sym(1, :), ...
                               coef_f_q_sym(2, :), ...
                               coef_f_q_sym(3, :), ...
                               coef_f_q_sym(4, :), ...
                               q));
t1_func = t_funcs{1};
t2_func = t_funcs{2};
t3_func = t_funcs{3};
mu_ = mu;

qq0 = q0;
last_q = 0;
last_err = 1e15;
for j = 1 : max_iter
    residual = abs(last_q) - abs(qq0);
    err = sqrt(residual.' * residual);
    if(err < 1e-20)
        break;
    elseif(err >= 1e-20 && last_err < err)
        mu_ = mu_ / 3;
    end
    last_err = err;
    
    last_q = qq0;
    Jaco = Jacob_fun(qq0);
    res = eq_fun(qq0);
    grad = Jaco.' * res;
    tmp = Jaco.' * Jaco;
    Hess = (tmp + mu_ * eye(4));
    qq0 = qq0 - Hess \ grad;
    qq0 = qq0 ./ norm(qq0);
end
    
t1_ = t1_func(pinvG, coefs_tq, qq0);
t2_ = t2_func(pinvG, coefs_tq, qq0);
t3_ = t3_func(pinvG, coefs_tq, qq0);
t = [t1_; t2_; t3_];

R = quat2dcm(qq0.');
X = [
    R, t;
    zeros(1, 3), 1];
end