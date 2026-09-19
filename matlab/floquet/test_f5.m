function test_f5()
% Test rapide : seulement f=5 (gamma=5), S=0.9, k de 0.1 a 40 par pas 0.1
% Ecrit S_0.90/f_5.0.csv et S_0.90/critical_test_f5.txt
% Utilise sigma_B_floquet et chebdif comme fonctions locales de
% B_freq_for_one_S.m -> on appelle donc la fonction principale avec un
% f_vec restreint via une copie locale du calcul.

addpath(pwd);

S    = 0.9;
eps  = 0.14;
E    = 0.05;
f    = 5.0;
om   = 2*f^2;
k_vec = 0.1:0.1:40.0;

out_dir = sprintf('S_%.2f', S);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

fprintf('Test f=%.1f, S=%.2f, om=%.1f\n', f, S, om);
fprintf('k : %g a %g (n=%d)\n\n', k_vec(1), k_vec(end), length(k_vec));

Wi_c_vec = nan(size(k_vec));
Wi_guess = 5;

t0 = tic;
for ik = 1:length(k_vec)
    k = k_vec(ik);
    try
        Wi_c = fzero(@(Wi) sigma_B_floquet(Wi, k, om, E, S, eps), Wi_guess);
        Wi_c = abs(Wi_c);                 % force la branche physique (positive)
        if Wi_c < 1e-3 || Wi_c > 200
            Wi_c_vec(ik) = NaN;
        else
            Wi_c_vec(ik) = Wi_c;
            Wi_guess     = Wi_c;
        end
    catch
        Wi_c_vec(ik) = NaN;
    end
    if mod(ik, 20) == 0
        fprintf('  ik=%d/%d  k=%.1f  Wi_c=%.3f  (%.1fs)\n', ...
            ik, length(k_vec), k, Wi_c_vec(ik), toc(t0));
    end
end
fprintf('Total: %.1fs\n', toc(t0));

% CSV
fname = fullfile(out_dir, sprintf('f_%.1f.csv', f));
fid = fopen(fname, 'w');
fprintf(fid, 'k, Wi_c\n');
for ik = 1:length(k_vec)
    fprintf(fid, '%g, %g\n', k_vec(ik), Wi_c_vec(ik));
end
fclose(fid);
fprintf('CSV ecrit : %s\n', fname);

% Critique
mask = Wi_c_vec > 0 & ~isnan(Wi_c_vec);
if any(mask)
    [Wic_min, j] = min(Wi_c_vec(mask));
    idxs = find(mask);
    k_c = k_vec(idxs(j));
    fprintf('Wi_c = %.4f  a  k_c = %.2f\n', Wic_min, k_c);
    fid = fopen(fullfile(out_dir, 'critical_test_f5.txt'), 'w');
    fprintf(fid, 'f\tWi_c\tk_c\n');
    fprintf(fid, '%g\t%g\t%g\n', f, Wic_min, k_c);
    fclose(fid);
end
end
