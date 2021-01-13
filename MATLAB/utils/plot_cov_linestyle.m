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


function plot_cov_linestyle(cov1, line1, cov2, line2, Data, mean_Data, point_style, i, j, circle_size, point_size)
global linestyle
AAA = cov1;
AA = [
    AAA(i, i), AAA(i, j);
    AAA(j, i), AAA(j, j);
    ];
linestyle = line1;
error_ellipse_linestyle(AA, 'scale', circle_size); hold on
plot(Data(i, :) - mean_Data(i), Data(j, :) - mean_Data(j), 'LineStyle', 'none', 'Marker', point_style, 'MarkerSize', point_size, 'Color', [0.5 0.5 0.5]); hold on
AAA = cov2;
AA = [
    AAA(i, i), AAA(i, j);
    AAA(j, i), AAA(j, j);
    ];
linestyle = line2;
error_ellipse_linestyle(AA, 'scale', circle_size); hold off
end
