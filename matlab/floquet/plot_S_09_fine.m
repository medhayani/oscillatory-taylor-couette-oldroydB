function plot_S_0_9_fine()
% Trace le resultat fin S=0.9 a partir du CSV deja produit

data = readmatrix('critical_Wi_vs_gamma_S_0.90.csv', 'NumHeaderLines', 1);
gamma = data(:,1); Wic = data(:,2); kc = data(:,3);

[Wic_min, imin] = min(Wic);
gamma_min = gamma(imin); kc_min = kc(imin);

figure('Position', [100 100 900 600]);
yyaxis left
plot(gamma, Wic, 'b-', 'LineWidth', 2); hold on;
plot(gamma_min, Wic_min, 'k*', 'MarkerSize', 16, 'LineWidth', 2);
ylabel('Wi_c', 'FontSize', 13);
yyaxis right
plot(gamma, kc, 'r--', 'LineWidth', 1.5);
ylabel('k_c', 'FontSize', 13);
xlabel('\gamma', 'FontSize', 13);
title(sprintf('S=0.9 grille fine : min Wi_c=%.3f a \\gamma=%.2f, k=%.2f', Wic_min, gamma_min, kc_min), 'FontSize', 13);
grid on;
saveas(gcf, 'fine_S_0.9_Wic_kc_vs_gamma.png');
saveas(gcf, 'fine_S_0.9_Wic_kc_vs_gamma.fig');

fprintf('Min Wi_c = %.4f a gamma = %.2f, k_c = %.2f\n', Wic_min, gamma_min, kc_min);
end
