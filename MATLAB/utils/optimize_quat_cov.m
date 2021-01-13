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



function [XX, ff] = optimize_quat_cov(AA, bb, q, F, cov_left, epsX, epsQuat, solver, verbose)
X = sdpvar(4, 4);
res = AA * vec(X) - bb;
f = res.' * res;
cons = [
    X - epsX * eye(4) >= 0, ...
    vec(F * X * F.' - cov_left) == 0, ...
    q.' * X * q <= epsQuat
    ];
options = sdpsettings('solver', solver, 'verbose', verbose);
if(strcmp(solver, 'sdpa_gmp'))
    options.sdpa_gmp.epsilonDash = 1.0e-35;
    options.sdpa_gmp.precision = 250;
end
optimize(cons, f, options);
XX = value(X);
ff = value(f);
end