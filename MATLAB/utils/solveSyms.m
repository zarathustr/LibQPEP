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



function sols = solveSyms(syst, vars)
A = jacobian(syst, vars);
rem = expand(syst - A * vars);
if(verLessThan('matlab', '8.0.0'))
    pinvA = inv(A);
else
    assumeAlso(A, 'real');
    pinvA = pinv(A);
end
sols_ = - pinvA * rem;

for i = 1 : length(sols_)
    str = sprintf('sols.%s = sols_(%d);', char(vars(i)), i);
    eval(str);
end
end