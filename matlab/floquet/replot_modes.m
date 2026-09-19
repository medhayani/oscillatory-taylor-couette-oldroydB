function replot_modes()
% Relit mode_characterization.csv et trace la figure avec legende correcte
% (harmonic / subharmonic / quasi-periodic), sans recalculer.

clear all; close all;

T = readtable('mode_characterization.csv');
gam = T.gamma; Wic = T.Wi_c; arg_m = T.arg_mu; sv_r = T.stress_vel_ratio; xmax = T.x_max;

% Classification rigoureuse base sur arg(mu)
harm = abs(arg_m) < 0.15;
sub  = abs(abs(arg_m) - pi) < 0.15;
qp   = ~harm & ~sub;

S = 0.9;

figure('Position', [100 100 1400 800]);

% --- Plot 1 : Wi_c colore par type ---
subplot(2,2,1);
plot(gam, Wic, 'k-', 'LineWidth', 1, 'HandleVisibility','off'); hold on;
h1 = plot(gam(harm), Wic(harm), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
h2 = plot(gam(sub),  Wic(sub),  'ms', 'MarkerSize', 7, 'MarkerFaceColor', 'm');
h3 = plot(gam(qp),   Wic(qp),   'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
xlabel('\gamma', 'FontSize', 12); ylabel('Wi_c', 'FontSize', 12);
% Construire legende uniquement avec les types presents
hh = [];
ll = {};
if any(harm), hh(end+1) = h1; ll{end+1} = sprintf('harmonique (%d pts)', sum(harm)); end
if any(sub),  hh(end+1) = h2; ll{end+1} = sprintf('sous-harmonique (%d pts)', sum(sub)); end
if any(qp),   hh(end+1) = h3; ll{end+1} = sprintf('quasi-periodique (%d pts)', sum(qp)); end
legend(hh, ll, 'Location', 'best');
title('Wi_c et type Floquet');
grid on;

% --- Plot 2 : arg(mu)/pi ---
subplot(2,2,2);
plot(gam(harm), arg_m(harm)/pi, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g'); hold on;
plot(gam(sub),  arg_m(sub)/pi,  'ms', 'MarkerSize', 7, 'MarkerFaceColor', 'm');
plot(gam(qp),   arg_m(qp)/pi,   'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
yline(0, 'g--', 'harmonique');
yline(1, 'm--', 'sub-harm.');
yline(-1, 'm--', 'sub-harm.');
yline(0.5, 'r:', 'arg=\pi/2');
yline(-0.5, 'r:', 'arg=-\pi/2');
xlabel('\gamma'); ylabel('arg(\mu)/\pi');
title('Argument du multiplicateur de Floquet');
grid on; ylim([-1.05, 1.05]);

% --- Plot 3 : ratio stress/velocity ---
subplot(2,2,3);
semilogy(gam(harm), sv_r(harm), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g'); hold on;
semilogy(gam(sub),  sv_r(sub),  'ms', 'MarkerSize', 7, 'MarkerFaceColor', 'm');
semilogy(gam(qp),   sv_r(qp),   'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
yline(1, 'k--');
xlabel('\gamma'); ylabel('||\tau|| / ||u||');
title('Ratio stress/vitesse dans le vecteur propre');
grid on;

% --- Plot 4 : localisation x_max ---
subplot(2,2,4);
plot(gam(harm), xmax(harm), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g'); hold on;
plot(gam(sub),  xmax(sub),  'ms', 'MarkerSize', 7, 'MarkerFaceColor', 'm');
plot(gam(qp),   xmax(qp),   'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
yline(0.5, 'k--', 'centre gap');
yline(0.2, 'k:', 'paroi int.');
yline(0.8, 'k:', 'paroi ext.');
xlabel('\gamma'); ylabel('x_{max} (max(|U|))');
ylim([0 1]);
title('Localisation spatiale du mode');
grid on;

sgtitle(sprintf('Caracterisation des modes critiques, S = %.1f (legende corrigee)', S), 'FontSize', 14);

saveas(gcf, 'mode_characterization_summary.png');
saveas(gcf, 'mode_characterization_summary.fig');

% Resume console
fprintf('\nClassification :\n');
fprintf('  Harmoniques     (arg=0)       : %d points\n', sum(harm));
fprintf('  Sous-harmoniques (arg=+/-pi)   : %d points\n', sum(sub));
fprintf('  Quasi-periodiques (autres)    : %d points\n', sum(qp));
fprintf('\nFigure regeneree : mode_characterization_summary.png\n');
end
