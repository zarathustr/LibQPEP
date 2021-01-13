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



function [ProjectedPoints, s] = generateProjectedPoints_(world_pt, K, R, t)
numPoints = size(world_pt, 1);
ProjectedPoints = zeros(numPoints, 2);
s = zeros(numPoints);
for i = 1 : numPoints
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates

    cameraMatrix = [R; t'] * K;
    projectedPoint = world_point * cameraMatrix;
    s(i) = projectedPoint(3);
    projectedPoint = projectedPoint(1 : 2) ./ s(i);
    ProjectedPoints(i, :) = projectedPoint;
end
end