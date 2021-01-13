function plot_R_cov(Sigma_R_stat, cov_R, Rs)
Lies = dcm2Lie(Rs).';
mean_Lie = mean(Lies, 2);
colors = linspecer(3);
subplot(2, 2, 1);
plot_cov_linestyle(Sigma_R_stat, {'--', 3, colors(1, :)}, cov_R, {'-', 2, colors(3, :)}, Lies, mean_Lie, '.', 1, 2, 3, 10.5);
xlabel('$\xi_1$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$\xi_2$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(2, 2, 2);
plot_cov_linestyle(Sigma_R_stat, {'--', 3, colors(1, :)}, cov_R, {'-', 2, colors(3, :)}, Lies, mean_Lie, '.', 1, 3, 3, 10.5);
xlabel('$\xi_1$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$\xi_3$', 'Interpreter', 'LaTeX', 'FontSize', 15);

subplot(2, 2, 4);
plot_cov_linestyle(Sigma_R_stat, {'--', 3, colors(1, :)}, cov_R, {'-', 2, colors(3, :)}, Lies, mean_Lie, '.', 2, 3, 3, 10.5);
xlabel('$\xi_2$', 'Interpreter', 'LaTeX', 'FontSize', 15);
ylabel('$\xi_3$', 'Interpreter', 'LaTeX', 'FontSize', 15);

legend({'Stat Covariance', 'Data Points', 'Nguyen et al. Covariance'}, 'Interpreter', 'LaTeX', 'FontSize', 15);

Sigma_Lie_stat = covx(Lies.', Lies.')
end