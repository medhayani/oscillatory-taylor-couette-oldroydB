function B_for_one_S(S)
% =========================================================================
% B_FOR_ONE_S : pour un S fixe par l'utilisateur, calcul Wi_c(k) pour
%               chaque gamma puis extraction des valeurs critiques
%
% USAGE :
%   B_for_one_S()        % utilise S = 1.0 par defaut
%   B_for_one_S(0.5)     % calcul avec S = 0.5
%
% SORTIES (CSV) :
%   1) Pour chaque gamma : un fichier "Wi_vs_k_gamma_X.X_S_X.XX.csv"
%      contenant deux colonnes (k, Wi_c(k)).
%   2) A la fin : un fichier "critical_Wi_vs_gamma_S_X.XX.csv"
%      contenant trois colonnes (gamma, Wi_c_min, k_c).
%
% Parametres physiques fixes :
%   eps = 0.14 (gap ratio)
%   E   = 0.05 (nombre elastique interne)
%   Liste de gamma : [1, 2, 3, 4, 5, 6, 7, 8] (a editer si besoin)
%   Liste de k     : 2 a 12 par pas de 0.5 (a editer si besoin)
% =========================================================================

if nargin < 1
    S = 1.0;
end

format short eng;

% ===== Parametres fixes =====
eps        = 0.14;
E          = 0.05;
gamma_list = 2:0.1:7;     % pas 0.1, 51 valeurs
k_vec      = 2:0.1:12;    % pas 0.1, 101 valeurs

fprintf('\n========== B_for_one_S : S = %.2f ==========\n', S);
fprintf('eps = %.2f, E = %.4f\n', eps, E);
fprintf('gamma_list = [%s]\n', num2str(gamma_list));
fprintf('k_vec : %g a %g par pas de %g (n=%d)\n\n', ...
        k_vec(1), k_vec(end), k_vec(2)-k_vec(1), length(k_vec));

% Stockage des criticalites
critical_data = zeros(length(gamma_list), 3);    % colonnes : gamma, Wi_c, k_c

t_total = tic;

for ig = 1:length(gamma_list)
    gamma = gamma_list(ig);
    om    = 2*gamma^2;
    De    = E*om;
    fprintf('--- gamma = %.1f, om = %.2f, De = %.3f ---\n', gamma, om, De);

    % Boucle sur k avec continuation
    Wi_c_vec = zeros(size(k_vec));
    Wi_guess = 5;
    for ik = 1:length(k_vec)
        k = k_vec(ik);
        try
            Wi_c = fzero(@(Wi) sigma_B_floquet(Wi, k, om, E, S, eps), Wi_guess);
            Wi_c_vec(ik) = Wi_c;
            Wi_guess = Wi_c;
        catch
            Wi_c_vec(ik) = NaN;
        end
    end

    % CSV pour ce gamma : Wi_c en fonction de k
    fname = sprintf('Wi_vs_k_gamma_%.1f_S_%.2f.csv', gamma, S);
    fid = fopen(fname, 'w');
    fprintf(fid, 'k, Wi_c\n');
    for ik = 1:length(k_vec)
        fprintf(fid, '%g, %g\n', k_vec(ik), Wi_c_vec(ik));
    end
    fclose(fid);

    % Extraction critique (robuste : gere le cas ou rien n'a converge)
    valid_mask = Wi_c_vec > 0 & ~isnan(Wi_c_vec);
    if any(valid_mask)
        valid = Wi_c_vec(valid_mask);
        [Wic_min, idx_in_valid] = min(valid);
        valid_indices = find(valid_mask);
        full_idx = valid_indices(idx_in_valid);
        k_c = k_vec(full_idx);
        critical_data(ig, :) = [gamma, Wic_min, k_c];
        fprintf('  -> Wi_c_min = %.4f a k_c = %.2f\n', Wic_min, k_c);
    else
        critical_data(ig, :) = [gamma, NaN, NaN];
        fprintf('  -> AUCUN Wi_c valide trouve (tous fzero ont echoue)\n');
    end
end

% CSV final consolide : critique en fonction de gamma
fname_crit = sprintf('critical_Wi_vs_gamma_S_%.2f.csv', S);
fid = fopen(fname_crit, 'w');
fprintf(fid, 'gamma, Wi_c, k_c\n');
for ig = 1:length(gamma_list)
    fprintf(fid, '%g, %g, %g\n', critical_data(ig,1), critical_data(ig,2), critical_data(ig,3));
end
fclose(fid);

fprintf('\n=== Termine en %.1f s ===\n', toc(t_total));
fprintf('Fichiers crees :\n');
fprintf('  - %d x Wi_vs_k_gamma_*_S_%.2f.csv\n', length(gamma_list), S);
fprintf('  - %s\n', fname_crit);

end
