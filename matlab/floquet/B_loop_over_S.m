function B_loop_over_S()
% =========================================================================
% B_LOOP_OVER_S : appelle B_for_one_S pour chaque valeur S de S_list,
%                 puis aggrege les fichiers de criticalite en un seul
%                 master CSV et trace Wi_c(gamma) pour tous les S.
%
% Pour chaque S :
%   - Appel a B_for_one_S(S) qui produit les CSV intermediaires.
%   - Lecture du fichier critical_Wi_vs_gamma_S_X.XX.csv produit.
%   - Concatenation dans la matrice all_critical.
%
% Sorties finales :
%   - master_critical_Wi.csv : colonnes (S, gamma, Wi_c, k_c)
%   - critical_Wi_vs_gamma_all_S.png/.fig : figure consolidee
% =========================================================================

clear all; close all; format short eng;

% ===== Liste de S a balayer (modifier si besoin) =====
S_list = [0.3, 0.5, 0.7, 0.9];

fprintf('\n========== B_loop_over_S : balayage de S ==========\n');
fprintf('S_list = [%s]\n\n', num2str(S_list));

t_total = tic;
all_critical = [];   % colonnes : S, gamma, Wi_c, k_c

for is = 1:length(S_list)
    S = S_list(is);
    fprintf('\n********** Lancement pour S = %.2f **********\n', S);
    B_for_one_S(S);

    % Lecture du fichier de criticalite produit
    fname_crit = sprintf('critical_Wi_vs_gamma_S_%.2f.csv', S);
    data = readmatrix(fname_crit, 'NumHeaderLines', 1);   % colonnes gamma, Wi_c, k_c
    n = size(data, 1);
    S_col = S * ones(n, 1);
    all_critical = [all_critical; [S_col, data]];
end

% ===== Master CSV =====
fid = fopen('master_critical_Wi.csv', 'w');
fprintf(fid, 'S, gamma, Wi_c, k_c\n');
for i = 1:size(all_critical, 1)
    fprintf(fid, '%g, %g, %g, %g\n', all_critical(i,1), all_critical(i,2), all_critical(i,3), all_critical(i,4));
end
fclose(fid);

% ===== Figure : Wi_c(gamma) pour chaque S =====
figure('Position', [100 100 900 600]);
colors = lines(length(S_list));
for is = 1:length(S_list)
    S = S_list(is);
    rows = all_critical(:,1) == S;
    gamma_S = all_critical(rows, 2);
    Wi_c_S  = all_critical(rows, 3);
    plot(gamma_S, Wi_c_S, 'o-', 'Color', colors(is,:), 'LineWidth', 1.8, ...
         'MarkerFaceColor', colors(is,:), 'MarkerSize', 8, ...
         'DisplayName', sprintf('S = %.2f', S));
    hold on;
end
xlabel('\gamma (frequence number)', 'FontSize', 13);
ylabel('Wi_c (critical Weissenberg)', 'FontSize', 13);
title('Systeme B : Wi_c critique vs \gamma pour plusieurs S', 'FontSize', 13);
legend('Location', 'best'); grid on;
saveas(gcf, 'critical_Wi_vs_gamma_all_S.png');
saveas(gcf, 'critical_Wi_vs_gamma_all_S.fig');

% Figure complementaire : k_c(gamma) pour chaque S
figure('Position', [100 100 900 600]);
for is = 1:length(S_list)
    S = S_list(is);
    rows = all_critical(:,1) == S;
    gamma_S = all_critical(rows, 2);
    k_c_S   = all_critical(rows, 4);
    plot(gamma_S, k_c_S, 's-', 'Color', colors(is,:), 'LineWidth', 1.8, ...
         'MarkerFaceColor', colors(is,:), 'MarkerSize', 8, ...
         'DisplayName', sprintf('S = %.2f', S));
    hold on;
end
xlabel('\gamma', 'FontSize', 13);
ylabel('k_c (critical wavenumber)', 'FontSize', 13);
title('Systeme B : k_c critique vs \gamma pour plusieurs S', 'FontSize', 13);
legend('Location', 'best'); grid on;
saveas(gcf, 'critical_kc_vs_gamma_all_S.png');
saveas(gcf, 'critical_kc_vs_gamma_all_S.fig');

fprintf('\n=========================================\n');
fprintf('=== B_loop_over_S termine en %.1f s ===\n', toc(t_total));
fprintf('=========================================\n');
fprintf('Fichiers consolides :\n');
fprintf('  - master_critical_Wi.csv\n');
fprintf('  - critical_Wi_vs_gamma_all_S.png/.fig\n');
fprintf('  - critical_kc_vs_gamma_all_S.png/.fig\n');

end
