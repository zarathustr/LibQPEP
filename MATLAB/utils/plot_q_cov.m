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



function plot_q_cov(Sigma_q_stat, cov_q, qs)

mean_q = mean(qs, 2);
colors = linspecer(3);
subplot(3, 3, 1);
plot_cov_linestyle(Sigma_q_stat, {'--', 3, colors(1, :)}, cov_q, {'-', 2, colors(2, :)}, qs, mean_q, '.', 1, 2, 3, 5.5);
xlabel('$q_0$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$q_1$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(3, 3, 2);
plot_cov_linestyle(Sigma_q_stat, {'--', 3, colors(1, :)}, cov_q, {'-', 2, colors(2, :)}, qs, mean_q, '.', 1, 3, 3, 5.5);
xlabel('$q_0$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$q_2$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(3, 3, 3);
plot_cov_linestyle(Sigma_q_stat, {'--', 3, colors(1, :)}, cov_q, {'-', 2, colors(2, :)}, qs, mean_q, '.', 1, 4, 3, 5.5);
xlabel('$q_0$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$q_3$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(3, 3, 5);
plot_cov_linestyle(Sigma_q_stat, {'--', 3, colors(1, :)}, cov_q, {'-', 2, colors(2, :)}, qs, mean_q, '.', 2, 3, 3, 5.5);
xlabel('$q_1$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$q_2$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(3, 3, 6);
plot_cov_linestyle(Sigma_q_stat, {'--', 3, colors(1, :)}, cov_q, {'-', 2, colors(2, :)}, qs, mean_q, '.', 2, 4, 3, 5.5);
xlabel('$q_1$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$q_3$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(3, 3, 9);
plot_cov_linestyle(Sigma_q_stat, {'--', 3, colors(1, :)}, cov_q, {'-', 2, colors(2, :)}, qs, mean_q, '.', 3, 4, 3, 5.5);
xlabel('$q_2$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$q_3$', 'Interpreter', 'LaTeX', 'FontSize', 15);

legend({'Stat Covariance of Quaternion', 'Data Points', 'Covariance from QPEP Optimization'}, 'Interpreter', 'LaTeX', 'FontSize', 14);
end