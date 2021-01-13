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



function [XX, ff] = optimize_quat_cov2(q, F, cov_left, solver, verbose)
if(~strcmp(solver, 'scs'))
    X = sdpvar(4, 4);
    f = q.' * X * q;
    cons = [
        X >= 0, ...
        sym2vec3(F * X * F.') == sym2vec3(cov_left)
        ];
    options = sdpsettings('solver', solver, 'verbose', verbose);
    if(strcmp(solver, 'sdpa_gmp'))
        options.sdpa_gmp.epsilonDash = 1.0e-35;
        options.sdpa_gmp.precision = 250;
    end
    optimize(cons, f, options);
    XX = value(X);
    ff = value(f);
else
    cvx_solver scs
    cvx_precision high
    
    cvx_begin sdp
         variable X(4,4) symmetric semidefinite
         dual variable Q
         minimize(q.' * X * q)
         X >= 0 * eye(4) : Q
         subject to
             sym2vec3(F * X * F.') == sym2vec3(cov_left)
    cvx_end
    XX = real(X);
    ff = 0;
end
end