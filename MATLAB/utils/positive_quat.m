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
% Select quaternions with positive scalar parts
% References:
%       [1]. Bar-Itzhack, I. Y. (2000). New Method for Extracting the Quaternion from a Rotation Matrix. 
%                      Journal of Guidance, Control, and Dynamics, 23(6), 1085–1087. 
%                      https://doi.org/10.2514/2.4654
%
%       [2]. Wu, J. (2019). Optimal Continuous Unit Quaternions from Rotation Matrices. 
%                      Journal of Guidance, Control, and Dynamics, 42(4), 919–922. 
%                      https://doi.org/10.2514/1.g004043


function qq = positive_quat(q)
if(q(1) < 0)
     qq = - q;
else
     qq = q;
end
end