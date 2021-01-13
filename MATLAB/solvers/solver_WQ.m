% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs), 
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu,  J.,  Zheng,  Y.,  Gao,  Z.,  Jiang,  Y.,  Hu,  X.,  Zhu,  Y.,  Jiao,  J.,  Liu,  M. (2020)
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


function sols = solver_WQ(data)
[C0, C1] = setup_elimination_template_solver_WQ(data);
C1 = C0 \ C1;
RR = [- C1(end - 31 : end, :); eye(40)];
AM_ind = [34, 1, 2, 3, 4, 5, 6, 7, 8, 35, 9, 10, 46, 11, 12, 13, 14, 15, 16, 17, 18, 47, 19, 20, 50, 21, 22, 23, 24, 25, 26, 51, 27, 28, 52, 29, 30, 53, 31, 32];
AM = RR(AM_ind, :);
[V, D] = eig(AM);
scale = sqrt(diag(D).' ./ ((V(13, :) .* V(35, :))));
V = V .* (ones(size(V, 1), 1) * scale);
sols = complex(zeros(5, 40));
sols(1, :) = V(1, :);
sols(2, :) = V(13, :);
sols(3, :) = V(25, :);
sols(4, :) = V(35, :);
sols(5, :) = V(3, :) ./ (sols(1, :) .* sols(2, :) .* sols(4, :));
end
