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
%
% Orthonormalization of matrices to SO(n)
% References: 
%      ﻿[1]. Arun, K. S., Huang, T. S., & Blostein, S. D. (1987). Least-Squares Fitting of Two 3-D Point Sets. 
%                IEEE Transactions on Pattern Analysis and Machine Intelligence, PAMI-9(5), 698–700.
%
%       [2]. ﻿Wu, J., Liu, M., Zhou, Z., & Li, R. (2020). Fast Symbolic 3-D Registration Solution. 
%                IEEE Transactions on Automation Science and Engineering, 17(2), 761–770. 
%                https://doi.org/10.1109/TASE.2019.2942324


function R = orthonormalize(A)
ss = size(A);
dim = ss(1);
[u, ~, v] = svd(A);
s = eye(dim, dim);
s(dim, dim) = det(u * v);
R = u * s * v';
end