function B_freq_for_one_SV1(S)
% =========================================================================
% B_FREQ_FOR_ONE_S : pour un S fixe, balaye la frequence et k, calcule
%                    Wi_c(k) puis extrait les valeurs critiques.
%                    FICHIER AUTONOME (sigma_B_floquet et chebdif inclus).
%
% USAGE :
%   B_freq_for_one_S()       % S = 1.0 par defaut
%   B_freq_for_one_S(0.7)    % calcul avec S = 0.7
%
% ORGANISATION DES SORTIES :
%   - Cree un dossier nomme "S_X.XX" dans le repertoire courant.
%   - Pour chaque frequence f : un CSV "f_X.X.csv"
%       colonnes (k, Wi_c) : seuil neutre principal (premier zero de sigma).
%   - A la fin : un TXT "critical_Wi_vs_f_S_X.XX.txt"
%       colonnes (f, Wi_c, k_c) ou Wi_c = min sur k du seuil neutre.
%
% Parametres balayes :
%   f : 0.5 -> 8.0 par pas de 0.1
%   k : 0.1 -> 40.0 par pas de 0.1
%
% Parametres physiques fixes :
%   eps = 0.14
%   E   = 0.05
% =========================================================================

if nargin < 1
    S = 0.3;
end

format short eng;

% ===== Parametres fixes =====
eps   = 0.14;
E     = 1;
f_vec = 0.5:0.1:8.0;     % frequence
k_vec = 0.1:0.1:40.0;    % k

Wi_guess0 = 5;           % devine initiale de fzero (continuation ensuite)

% ===== Dossier de sortie =====
out_dir = sprintf('S_%.2f', S);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('\n========== B_freq_for_one_S : S = %.2f ==========\n', S);
fprintf('eps = %.2f, E = %.4f\n', eps, E);
fprintf('f : %g a %g par pas de %g (n=%d)\n', ...
        f_vec(1), f_vec(end), f_vec(2)-f_vec(1), length(f_vec));
fprintf('k : %g a %g par pas de %g (n=%d)\n', ...
        k_vec(1), k_vec(end), k_vec(2)-k_vec(1), length(k_vec));
fprintf('Dossier de sortie : %s\n\n', out_dir);

critical_data = zeros(length(f_vec), 3);    % colonnes : f, Wi_c, k_c

t_total = tic;

for ifr = 1:length(f_vec)
    f  = f_vec(ifr);
    om = 2*f^2;        % meme relation que dans B_for_one_S (gamma -> om)
    De = E*om;
    fprintf('--- f = %.1f, om = %.2f, De = %.3f ---\n', f, om, De);

    % Boucle sur k avec continuation (un seul fzero par k -> rapide)
    Wi_c_vec = nan(size(k_vec));
    Wi_guess = Wi_guess0;
    for ik = 1:length(k_vec)
        k = k_vec(ik);
        try
            Wi_c = fzero(@(Wi) sigma_B_floquet(Wi, k, om, E, S, eps), Wi_guess);
            Wi_c = abs(Wi_c);                 % branche physique (positive)
            if Wi_c < 1e-3 || Wi_c > 200
                Wi_c_vec(ik) = NaN;
            else
                Wi_c_vec(ik) = Wi_c;
                Wi_guess     = Wi_c;
            end
        catch
            Wi_c_vec(ik) = NaN;
        end
    end

    % CSV pour cette frequence : (k, Wi_c)
    fname = fullfile(out_dir, sprintf('f_%.1f.csv', f));
    fid   = fopen(fname, 'w');
    fprintf(fid, 'k, Wi_c\n');
    for ik = 1:length(k_vec)
        fprintf(fid, '%g, %g\n', k_vec(ik), Wi_c_vec(ik));
    end
    fclose(fid);

    % Extraction critique : min sur k de Wi_c
    valid_mask = Wi_c_vec > 0 & ~isnan(Wi_c_vec);
    if any(valid_mask)
        valid                 = Wi_c_vec(valid_mask);
        [Wic_min, idx_in_val] = min(valid);
        valid_indices         = find(valid_mask);
        full_idx              = valid_indices(idx_in_val);
        k_c                   = k_vec(full_idx);
        critical_data(ifr, :) = [f, Wic_min, k_c];
        fprintf('  -> Wi_c_min = %.4f a k_c = %.2f\n', Wic_min, k_c);
    else
        critical_data(ifr, :) = [f, NaN, NaN];
        fprintf('  -> AUCUN Wi_c valide trouve\n');
    end
end

% TXT final consolide
fname_crit = fullfile(out_dir, sprintf('critical_Wi_vs_f_S_%.2f.txt', S));
fid = fopen(fname_crit, 'w');
fprintf(fid, 'f\tWi_c\tk_c\n');
for ifr = 1:length(f_vec)
    fprintf(fid, '%g\t%g\t%g\n', ...
            critical_data(ifr,1), critical_data(ifr,2), critical_data(ifr,3));
end
fclose(fid);

fprintf('\n=== Termine en %.1f s ===\n', toc(t_total));
fprintf('Fichiers crees dans %s :\n', out_dir);
fprintf('  - %d x f_*.csv\n', length(f_vec));
fprintf('  - %s\n', fname_crit);

end


% =========================================================================
% ============ FONCTIONS LOCALES (systeme B + Chebyshev) ==================
% =========================================================================

function sigma = sigma_B_floquet(Wi, k, om, E, S, eps)
% Plus grand exposant de Floquet (partie reelle) pour le systeme B
% (formulation visqueuse corrigee avec terme de courbure elastique)

T = 2*pi/om;
tol = 1e-2;
N = 8;
phi0 = eye(8*N-6);
x0 = reshape(phi0, [], 1);
options = odeset('RelTol', tol, 'AbsTol', tol);
[~, F] = ode45(@(t, F) rhs_B(t, F, Wi, k, om, E, S, eps, N), [0 T], x0, options);
Mon = reshape(F(end, :), 8*N-6, 8*N-6);
sig_eig = eig(Mon);
ee = log(sig_eig)/T;
[~, is] = sort(-real(ee));
sigma = real(ee(is(1)));
end


function Fdot = rhs_B(t, F, Wi, k, om, E, S, eps, N)
% RHS de d/dt phi = qqq * phi (apres elimination des CL)

[x, DM] = chebdif(N, 4);
s = 2;
D1 = s*DM(:,:,1); D2 = s^2*DM(:,:,2); D4 = s^4*DM(:,:,4);
i = sqrt(-1); I = eye(N); Z = zeros(N);

% Etat de base
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

% Matrice B (LHS)
B = zeros(8*N);
B(3:N-2,:) = [D2(3:N-2,1:N) - k^2*I(3:N-2,1:N), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:)];
B(N+2:2*N-1,:) = [Z(2:N-1,1:N), I(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N)];
for blk = 3:8
    B((blk-1)*N+1:blk*N, (blk-1)*N+1:blk*N) = E*I;
end

% Matrice A (RHS) avec courbure elastique
A = zeros(8*N);
A(3:N-2,:) = [(1-S)*(D4(3:N-2,:) - 2*k^2*D2(3:N-2,:) + k^4*I(3:N-2,:)), Z(3:N-2,:), -k^2*D1(3:N-2,:), Z(3:N-2,:), -i*k^3*I(3:N-2,:) - i*k*D2(3:N-2,:), eps*k^2*I(3:N-2,:), Z(3:N-2,:), k^2*D1(3:N-2,:)];
A(N+2:2*N-1,:) = [Z(2:N-1,1:N), (1-S)*(D2(2:N-1,:) - k^2*I(2:N-1,:)), Z(2:N-1,:), D1(2:N-1,:) + 2*eps*I(2:N-1,:), Z(2:N-1,:), Z(2:N-1,:), i*k*I(2:N-1,:), Z(2:N-1,:)];
A(2*N+1:3*N,:) = [2*S*D1, Z, -I, Z, Z, Z, Z, Z];
A(3*N+1:4*N,:) = [Wi*Tc*D1 - Wi*dTs*I + eps*Wi*Tc*I, S*D1, Wi*dVs*I - eps*Wi*Vc*I, -I, Z, Z, Z, Z];
A(4*N+1:5*N,:) = [i*(1/k)*S*D2 + i*k*S*I, Z, Z, Z, -I, Z, Z, Z];
A(5*N+1:6*N,:) = [-Wi*dTss*I + 2*eps*Wi*Tcc*I, 2*Wi*Tc*D1 - 2*eps*Wi*Tc*I, Z, 2*Wi*dVs*I - 2*eps*Wi*Vc*I, Z, -I, Z, Z];
A(6*N+1:7*N,:) = [i*(1/k)*Wi*(Tc*D2 + eps*Tc*D1 - eps^2*Tc*I), i*k*S*I, Z, Z, Wi*dVs*I - eps*Wi*Vc*I, Z, -I, Z];
A(7*N+1:8*N,:) = [-2*S*D1 - 2*eps*S*I, Z, Z, Z, Z, Z, Z, -I];

% Conditions aux limites
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


function [x, DM] = chebdif(N, M)
% Matrices de differentiation Chebyshev-Gauss-Lobatto (Weideman-Reddy)

I = eye(N); L = logical(I);
n1 = floor(N/2); n2 = ceil(N/2);
k = (0:N-1)'; th = k*pi/(N-1);
x = sin(pi*(N-1:-2:1-N)'/(2*(N-1)));
T = repmat(th/2, 1, N);
DX = 2*sin(T'+T).*sin(T'-T);
DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))]; DX(L) = ones(N,1);
C = toeplitz((-1).^k);
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2; C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;
Z = 1./DX; Z(L) = zeros(N,1);
D = eye(N);
for ell = 1:M
    D = ell*Z.*(C.*repmat(diag(D),1,N) - D);
    D(L) = -sum(D');
    DM(:,:,ell) = D;
end
end
