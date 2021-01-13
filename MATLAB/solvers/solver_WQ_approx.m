% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),   
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu,    J.,    Zheng,    Y.,    Gao,    Z.,    Jiang,    Y.,    Hu,    X.,    Zhu,    Y.,    Jiao,    J.,    Liu,    M. (2020)
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


function sols = solver_WQ_approx(data)
[C0, C1_] = setup_elimination_template_solver_WQ_approx(data);
C1 = C0 \ C1_;
RR = [ - C1(end-12:end, :); eye(27)];
AM_ind = [38, 16, 1, 19, 2, 3, 21, 22, 4, 25, 5, 6, 28, 7, 8, 30, 31, 9, 34, 10, 11, 36, 37, 12, 39, 40, 13];
AM = RR(AM_ind, :);
[V, D] = eig(AM);
V = V ./ (ones(size(V, 1), 1) * V(1, :));
sols = complex(zeros(4, 27));
sols(2, :) = sqrt(V(10, :));
sols(1, :) = V(2, :) ./ (sols(2, :));
sols(3, :) = V(4, :) ./ (sols(1, :));
sols(4, :) = V(3, :) ./ (sols(1, :) .* V(16, :));
end



