function characterize_modes()
% =========================================================================
% CHARACTERIZE_MODES : pour chaque point critique (gamma, Wi_c) du
% fichier critical_Wi_vs_gamma_S_0.90.csv (grille fine), calcule la
% caracterisation complete du mode :
%
%   - mu_c, |mu|, arg(mu)               : multiplicateur de Floquet critique
%   - type_floquet                      : harmonic / subharmonic / quasi-periodic
%   - stress_vel_ratio                  : ||tau||/||u|| dans le vecteur propre
%   - x_max                             : localisation spatiale du max de |U|
%   - dom_stress                        : composante de stress dominante
%   - phi_ela_max                       : max(|Phi_ela|) au point critique
%   - nature                            : interpretation physique
%
% Sortie : mode_characterization.csv
% =========================================================================

clear all; close all; format short eng;

% ===== Parametres =====
S      = 0.9;
eps    = 0.14;
E      = 0.05;
csv_in = 'critical_Wi_vs_gamma_S_0.90.csv';

% Lecture des points critiques
data = readmatrix(csv_in, 'NumHeaderLines', 1);
gamma_vec = data(:,1); Wi_c_vec = data(:,2); k_c_vec = data(:,3);
n_points = length(gamma_vec);

fprintf('\n========== CARACTERISATION DES MODES ==========\n');
fprintf('Lecture de %s : %d points critiques\n', csv_in, n_points);
fprintf('S = %.2f, eps = %.2f, E = %.4f\n\n', S, eps, E);

% Stockage des resultats
results = struct();

fprintf('%6s | %10s | %10s | %10s | %10s | %12s | %8s | %18s\n', ...
        'gamma', 'Wi_c', 'k_c', '|mu|', 'arg(mu)', 'type', 'tau/u', 'nature');
fprintf('-------|------------|------------|------------|------------|--------------|----------|--------------------\n');

t_start = tic;
for i = 1:n_points
    gamma = gamma_vec(i);
    Wi_c  = Wi_c_vec(i);
    k_c   = k_c_vec(i);
    if isnan(Wi_c) || isnan(k_c)
        fprintf('%6.2f | %10s\n', gamma, 'NaN');
        continue;
    end

    om = 2*gamma^2;
    [mu_c, type_fl, sv_ratio, x_max, dom_str, phi_ela_max, nature] = ...
        characterize_one(Wi_c, k_c, om, E, S, eps);

    results(i).gamma = gamma;
    results(i).Wi_c = Wi_c;
    results(i).k_c = k_c;
    results(i).mu_c = mu_c;
    results(i).abs_mu = abs(mu_c);
    results(i).arg_mu = angle(mu_c);
    results(i).type = type_fl;
    results(i).sv_ratio = sv_ratio;
    results(i).x_max = x_max;
    results(i).dom_str = dom_str;
    results(i).phi_ela_max = phi_ela_max;
    results(i).nature = nature;

    fprintf('%6.2f | %10.4f | %10.2f | %10.4f | %10.4f | %12s | %8.2f | %18s\n', ...
            gamma, Wi_c, k_c, abs(mu_c), angle(mu_c), type_fl, sv_ratio, nature);
end

% Sauvegarde CSV
fid = fopen('mode_characterization.csv', 'w');
fprintf(fid, 'gamma, Wi_c, k_c, abs_mu, arg_mu, type_floquet, stress_vel_ratio, x_max, dominant_stress, phi_ela_max, nature\n');
for i = 1:n_points
    if isfield(results(i), 'gamma') && ~isempty(results(i).gamma)
        r = results(i);
        fprintf(fid, '%g, %g, %g, %g, %g, %s, %g, %g, %s, %g, %s\n', ...
                r.gamma, r.Wi_c, r.k_c, r.abs_mu, r.arg_mu, r.type, r.sv_ratio, r.x_max, r.dom_str, r.phi_ela_max, r.nature);
    end
end
fclose(fid);

fprintf('\n=== Termine en %.1f s ===\n', toc(t_start));
fprintf('Fichier : mode_characterization.csv\n');

% Plot consolide
make_summary_plot(results, S);

end


function [mu_c, type_fl, sv_ratio, x_max, dom_str, phi_ela_max, nature] = characterize_one(Wi, k, om, E, S, eps)
% Pour un (Wi_c, k_c, gamma), calcule les caracteristiques du mode dominant

T = 2*pi/om;
tol = 1e-3;
N = 8;
n_dim = 8*N - 6;

phi0 = eye(n_dim);
x0 = reshape(phi0, [], 1);
options = odeset('RelTol', tol, 'AbsTol', tol);
[~, F] = ode45(@(t, F) rhs_B(t, F, Wi, k, om, E, S, eps, N), [0 T], x0, options);
Mon = reshape(F(end, :), n_dim, n_dim);

% Decomposition spectrale
[V, D] = eig(Mon);
mu_all = diag(D);

% Multiplicateur le plus proche du cercle unite
[~, idx] = min(abs(abs(mu_all) - 1));
mu_c = mu_all(idx);
eigvec = V(:, idx);

% Type Floquet
arg_mu = angle(mu_c);
if abs(arg_mu) < 0.15
    type_fl = 'harmonic';
elseif abs(abs(arg_mu) - pi) < 0.15
    type_fl = 'subharmonic';
else
    type_fl = 'quasi-period';
end

% Composantes du vecteur propre apres BC elimination :
% indices 1:N-4 = U(x_int), N-3:2N-6 = V(x_int), 2N-5:8N-6 = stresses (6 blocs de N)
n_vel = 2*N - 6;
U_int = abs(eigvec(1:N-4));
V_int = abs(eigvec(N-3:n_vel));
RR = abs(eigvec(n_vel+1:n_vel+N));
Rth = abs(eigvec(n_vel+N+1:n_vel+2*N));
RZ = abs(eigvec(n_vel+2*N+1:n_vel+3*N));
thth = abs(eigvec(n_vel+3*N+1:n_vel+4*N));
thZ = abs(eigvec(n_vel+4*N+1:n_vel+5*N));
ZZ = abs(eigvec(n_vel+5*N+1:n_vel+6*N));

norm_vel = sqrt(norm(U_int)^2 + norm(V_int)^2);
norm_str = sqrt(norm(RR)^2 + norm(Rth)^2 + norm(RZ)^2 + norm(thth)^2 + norm(thZ)^2 + norm(ZZ)^2);

if norm_vel < 1e-10
    sv_ratio = 1e10;
else
    sv_ratio = norm_str / norm_vel;
end

% Localisation spatiale du max de |U|
[~, x_idx] = max(U_int);
% U_int corresponds to x interior points 3 to N-2 (N-4 values)
% x_chebyshev grid : x = (1+xi)/2 with xi = cos(k pi/(N-1))
[xi_full, ~] = chebdif(N, 1);
x_full = 0.5*(1 + xi_full);
x_int = x_full(3:N-2);
x_max = x_int(x_idx);

% Composante de stress dominante
str_norms = [norm(RR), norm(Rth), norm(RZ), norm(thth), norm(thZ), norm(ZZ)];
str_names = {'RR', 'Rth', 'RZ', 'thth', 'thZ', 'ZZ'};
[~, dom_idx] = max(str_norms);
dom_str = str_names{dom_idx};

% Phi_ela = -E*S*dN1/dx maximum sur le cycle
[t_out, F_traj] = ode45(@(t, F) rhs_B(t, F, Wi, k, om, E, S, eps, N), [0 T], x0, options);
% On utilise la base evaluee a t = T/4 (max de cos pour Vi)
t_eval = T/4;
[V1, V2, T1, T2, T3, T4, T5, T6, T7, T8] = base_flow_components(om, E, S, N);
% N1 = thth - RR (mais RR est petit dans la base)
% On approxime N1 ~ thth
% dthth/dx ~ derivee spatiale du thth de base au temps T/4
% thth_base(x, t) = Wi*(eps*T7 + T8) + Wi*((eps*T3+T4)*cos(2omt) + (eps*T5+T6)*sin(2omt))
% Au temps t=0 par exemple :
thth_base = Wi*(eps*T7 + T8) + Wi*(eps*T3 + T4);
[~, DM] = chebdif(N, 1); D1 = 2*DM(:,:,1);
dthth = D1 * thth_base;
phi_ela_x = -E*S*dthth;
phi_ela_max = max(abs(phi_ela_x));

% Nature physique : heuristique base sur les indicateurs
% (en regime purement elastique Re=0, distinguer entre :
%  - localise paroi (x_max < 0.2 ou > 0.8) : couche de Stokes
%  - centre gap (0.4 < x_max < 0.6) : LSM bulk mode
%  - asymetrique : interaction couches )

if x_max < 0.2 || x_max > 0.8
    loc = 'Stokes-layer';
elseif x_max > 0.4 && x_max < 0.6
    loc = 'bulk-LSM';
else
    loc = 'asymmetric';
end

if strcmp(type_fl, 'subharmonic')
    par = 'param-resonance';
elseif strcmp(type_fl, 'harmonic')
    par = 'synchronous';
else
    par = 'incommensurate';
end

% Combinaison
nature = sprintf('%s-%s', loc, par);

end


function [V1, V2, T1, T2, T3, T4, T5, T6, T7, T8] = base_flow_components(om, E, S, N)
% Composantes spatiales de l'ecoulement de base

[x, DM] = chebdif(N, 1); D1 = 2*DM(:,:,1);
RR = (1-S)*E; gamma_v = sqrt(om/2);
C1 = (om*E)^2 + 1; C2 = (om*RR)^2 + 1; C3 = om*(E - RR);
beta = sqrt(sqrt(C1/C2) + C3/C2); zeta = sqrt(sqrt(C1/C2) - C3/C2);
xx = 0.5*(1+x); xxx = 0.5*(1-x); epp = 1;
C4 = epp*cos(beta*gamma_v) + cosh(zeta*gamma_v);
C7 = om*S*E; C8 = (2*om*E)^2 + 1;
V11 = cos(beta*gamma_v*xx).*cosh(zeta*gamma_v*xxx) + epp*cos(beta*gamma_v*xxx).*cosh(zeta*gamma_v*xx);
V22 = sin(beta*gamma_v*xx).*sinh(zeta*gamma_v*xxx) + epp*sin(beta*gamma_v*xxx).*sinh(zeta*gamma_v*xx);
V1 = V11./C4; V2 = V22./C4;
dV1 = D1*V1; dV2 = D1*V2;
T1 = (S/C1)*dV1 - (C7/C1)*dV2;
T2 = (C7/C1)*dV1 + (S/C1)*dV2;

% T3-T8 needs Wi but for base normal stress shape, use Wi=1 :
Wi_norm = 1;
C9 = Wi_norm; C10 = om*Wi_norm*E;
T3 = -(C9/C8)*T1.*V1 + 2*(C10/C8)*T1.*V2 + 2*(C10/C8)*T2.*V1 + (C9/C8)*T2.*V2;
T4 = (C9/C8)*T1.*dV1 - 2*(C10/C8)*T1.*dV2 - 2*(C10/C8)*T2.*dV1 - (C9/C8)*T2.*dV2;
T5 = -2*(C10/C8)*T1.*V1 - (C9/C8)*T1.*V2 - (C9/C8)*T2.*V1 - 2*(C10/C8)*T2.*V2;
T6 = 2*(C10/C8)*T1.*dV1 + (C9/C8)*T1.*dV2 + (C9/C8)*T2.*dV1 - 2*(C10/C8)*T2.*dV2;
T7 = -C9*T1.*V1 - C9*T2.*V2;
T8 = C9*T1.*dV1 + C9*T2.*dV2;
end


function make_summary_plot(results, S)
% Figure recapitulative : Wi_c et type Floquet vs gamma

valid_idx = arrayfun(@(r) ~isempty(r.gamma), results);
results = results(valid_idx);
gam = arrayfun(@(r) r.gamma, results);
Wic = arrayfun(@(r) r.Wi_c, results);
arg_m = arrayfun(@(r) r.arg_mu, results);
sv_r = arrayfun(@(r) r.sv_ratio, results);

figure('Position', [100 100 1200 700]);

subplot(2,2,1);
plot(gam, Wic, 'b-', 'LineWidth', 1.5); hold on;
% Coloriage par type Floquet
harm = abs(arg_m) < 0.15;
sub = abs(abs(arg_m) - pi) < 0.15;
plot(gam(harm), Wic(harm), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
plot(gam(sub), Wic(sub), 'ms', 'MarkerSize', 6, 'MarkerFaceColor', 'm');
plot(gam(~harm & ~sub), Wic(~harm & ~sub), 'r^', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
xlabel('\gamma'); ylabel('Wi_c');
legend('Wi_c(\gamma)', 'harmonic', 'subharmonic', 'quasi-periodic', 'Location', 'best');
title('Wi_c et type Floquet'); grid on;

subplot(2,2,2);
plot(gam, arg_m/pi, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('\gamma'); ylabel('arg(\mu)/\pi');
yline(0, 'g--'); yline(1, 'm--'); yline(-1, 'm--');
title('Argument du multiplicateur (en \pi)'); grid on;

subplot(2,2,3);
semilogy(gam, sv_r, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('\gamma'); ylabel('||\tau|| / ||u||');
title('Ratio stress/vitesse dans le vecteur propre'); grid on;

subplot(2,2,4);
xmax = arrayfun(@(r) r.x_max, results);
plot(gam, xmax, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('\gamma'); ylabel('x_{max}');
yline(0.5, 'k--'); ylim([0 1]);
title('Localisation spatiale du max(|U|)'); grid on;

sgtitle(sprintf('Caracterisation des modes critiques, S = %.1f', S), 'FontSize', 14);
saveas(gcf, 'mode_characterization_summary.png');
saveas(gcf, 'mode_characterization_summary.fig');
end


function Fdot = rhs_B(t, F, Wi, k, om, E, S, eps, N)
% Identique a sigma_B_floquet.m

[x, DM] = chebdif(N, 4);
s = 2;
D1 = s*DM(:,:,1); D2 = s^2*DM(:,:,2); D4 = s^4*DM(:,:,4);
i = sqrt(-1); I = eye(N); Z = zeros(N);

RR = (1-S)*E; gamma_v = sqrt(om/2);
C1 = (om*E)^2 + 1; C2 = (om*RR)^2 + 1; C3 = om*(E - RR);
beta = sqrt(sqrt(C1/C2) + C3/C2); zeta = sqrt(sqrt(C1/C2) - C3/C2);
xx = 0.5*(1+x); xxx = 0.5*(1-x); epp = 1;
C4 = epp*cos(beta*gamma_v) + cosh(zeta*gamma_v);
C7 = om*S*E; C8 = (2*om*E)^2 + 1;
V11 = cos(beta*gamma_v*xx).*cosh(zeta*gamma_v*xxx) + epp*cos(beta*gamma_v*xxx).*cosh(zeta*gamma_v*xx);
V22 = sin(beta*gamma_v*xx).*sinh(zeta*gamma_v*xxx) + epp*sin(beta*gamma_v*xxx).*sinh(zeta*gamma_v*xx);
V1 = V11./C4; V2 = V22./C4;
dV1 = D1*V1; dV2 = D1*V2;
T1 = (S/C1)*dV1 - (C7/C1)*dV2;
T2 = (C7/C1)*dV1 + (S/C1)*dV2;
C9_n = Wi; C10_n = om*Wi*E;
T3 = -(C9_n/C8)*T1.*V1 + 2*(C10_n/C8)*T1.*V2 + 2*(C10_n/C8)*T2.*V1 + (C9_n/C8)*T2.*V2;
T4 = (C9_n/C8)*T1.*dV1 - 2*(C10_n/C8)*T1.*dV2 - 2*(C10_n/C8)*T2.*dV1 - (C9_n/C8)*T2.*dV2;
T5 = -2*(C10_n/C8)*T1.*V1 - (C9_n/C8)*T1.*V2 - (C9_n/C8)*T2.*V1 - 2*(C10_n/C8)*T2.*V2;
T6 = 2*(C10_n/C8)*T1.*dV1 + (C9_n/C8)*T1.*dV2 + (C9_n/C8)*T2.*dV1 - 2*(C10_n/C8)*T2.*dV2;
T7 = -C9_n*T1.*V1 - C9_n*T2.*V2;
T8 = C9_n*T1.*dV1 + C9_n*T2.*dV2;

Vi = V1*cos(om*t) + V2*sin(om*t);
Ti = T1*cos(om*t) + T2*sin(om*t);
Tii = (eps*T3 + T4)*cos(2*om*t) + (eps*T5 + T6)*sin(2*om*t) + (eps*T7 + T8);
Vc = diag(Vi); dVs = diag(D1*Vi); Tc = diag(Ti); dTs = diag(D1*Ti);
Tcc = diag(Tii); dTss = diag(D1*Tii);

B = zeros(8*N);
B(3:N-2,:) = [D2(3:N-2,1:N) - k^2*I(3:N-2,1:N), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:)];
B(N+2:2*N-1,:) = [Z(2:N-1,1:N), I(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N)];
for blk = 3:8
    B((blk-1)*N+1:blk*N, (blk-1)*N+1:blk*N) = E*I;
end

A = zeros(8*N);
A(3:N-2,:) = [(1-S)*(D4(3:N-2,:) - 2*k^2*D2(3:N-2,:) + k^4*I(3:N-2,:)), Z(3:N-2,:), -k^2*D1(3:N-2,:), Z(3:N-2,:), -i*k^3*I(3:N-2,:) - i*k*D2(3:N-2,:), eps*k^2*I(3:N-2,:), Z(3:N-2,:), k^2*D1(3:N-2,:)];
A(N+2:2*N-1,:) = [Z(2:N-1,1:N), (1-S)*(D2(2:N-1,:) - k^2*I(2:N-1,:)), Z(2:N-1,:), D1(2:N-1,:) + 2*eps*I(2:N-1,:), Z(2:N-1,:), Z(2:N-1,:), i*k*I(2:N-1,:), Z(2:N-1,:)];
A(2*N+1:3*N,:) = [2*S*D1, Z, -I, Z, Z, Z, Z, Z];
A(3*N+1:4*N,:) = [Wi*Tc*D1 - Wi*dTs*I + eps*Wi*Tc*I, S*D1, Wi*dVs*I - eps*Wi*Vc*I, -I, Z, Z, Z, Z];
A(4*N+1:5*N,:) = [i*(1/k)*S*D2 + i*k*S*I, Z, Z, Z, -I, Z, Z, Z];
A(5*N+1:6*N,:) = [-Wi*dTss*I + 2*eps*Wi*Tcc*I, 2*Wi*Tc*D1 - 2*eps*Wi*Tc*I, Z, 2*Wi*dVs*I - 2*eps*Wi*Vc*I, Z, -I, Z, Z];
A(6*N+1:7*N,:) = [i*(1/k)*Wi*(Tc*D2 + eps*Tc*D1 - eps^2*Tc*I), i*k*S*I, Z, Z, Wi*dVs*I - eps*Wi*Vc*I, Z, -I, Z];
A(7*N+1:8*N,:) = [-2*S*D1 - 2*eps*S*I, Z, Z, Z, Z, Z, Z, -I];

c1 = [I(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:)];
c2 = [D1(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:)];
c3 = [I(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:)];
c4 = [D1(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:)];
c5 = [Z(1,:), I(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:)];
c6 = [Z(N,:), I(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:)];
C = [c1; c2; c3; c4; c5; c6];
r = [1, 2, N-1, N, N+1, 2*N];
ki = [3:N-2, N+2:2*N-1, 2*N+1:8*N];
G = -pinv(C(:,r))*C(:,ki);
AA = A(ki,ki) + A(ki,r)*G;
BB = B(ki,ki) + B(ki,r)*G;
qqq = AA/BB;

phi = reshape(F(1:end), 8*N-6, 8*N-6);
Fdot = reshape(qqq*phi, (8*N-6)^2, 1);
end
