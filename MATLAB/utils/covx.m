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


function cov = covx(X, Y)
len = size(X, 1);
mean_X = mean(X, 1);
mean_Y = mean(Y, 1);
nX = size(X, 2);
nY = size(Y, 2);
cov = zeros(nX, nY);
for i = 1 : len
    res = (X(i, :) - mean_X).' * (Y(i, :) - mean_Y);
    cov = cov + 1 / len * res;
end
end