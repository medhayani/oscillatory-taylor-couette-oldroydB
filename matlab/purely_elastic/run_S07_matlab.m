function run_S07_matlab(varargin)
% =========================================================================
% RUN_S07_MATLAB : Wi_c(f) a S = 0.7 (Oldroyd-B, terme visqueux (1-S) actif),
%                  systeme B / Floquet, N = 12, pour les 5 E de la campagne.
%                  FICHIER AUTONOME (moteur + chebdif inclus).
%
% Moteur ROBUSTE (identique au pipeline Python valide, PAS la version N=8/ode45) :
%   * matrices A/B = portage VERBATIM de rhs_B (sigma_B_floquet.m)
%   * monodromie par produit d'exponentielles (expm) avec RENORMALISATION a
%     chaque sous-pas -> pas d'overflow a Wi/E extremes (ou ode45 echoue)
%   * rejet du mode PARASITE longue-onde : si sigma(Wi->0) >= 0, le k est NaN
%   * Wi_c = 1er passage sigma:-/+ (seuil le plus bas) + fzero BRACKETE (pas de
%     saut de branche comme le fzero a devinette de B_freq_for_one_S)
%
% WI_HI RECALIBRE PAR E (a S=0.7 le seuil est ~x10 plus haut qu'a S=1 ;
%   sonde _calib_WIHI_S07). Sans ca -> tout NaN a E>=10.
%
% USAGE (sur le PC puissant) :
%   run_S07_matlab              % les 5 E : 0.006, 0.01, 1, 10, 50
%   run_S07_matlab(10, 50)      % seulement E=10 et E=50 (les lourds)
%   run_S07_matlab(1)           % seulement E=1
%
% PARALLELISME : parfor sur f. Avec Parallel Computing Toolbox -> multi-coeurs
%   (demarre un pool automatiquement). Sans PCT, parfor s'execute en serie.
%
% SORTIES (format identique au Python -> lisible par _ladder_final.py) :
%   ucm_E<tag>/S_0.70/ : f_<f>.csv (k, Wi_c)  +  critical_Wi_vs_f_S_0.70.txt
%   Reprise : un f_*.csv deja present est saute (mettre FORCE=true pour refaire).
%
% Auto-test au demarrage : verifie sigma vs 3 references du moteur Python
%   valide ; ABANDONNE si le portage est faux.
% =========================================================================

% -------------------- Parametres globaux --------------------
EPS   = 0.14;
Nc    = 12;                 % points Chebyshev
M_SUB = 40;                 % sous-pas expm / periode
S     = 0.70;
WI_LO = 0.05;
FORCE = false;              % true = recalcule meme si le f_*.csv existe

% ---- Grille de FREQUENCES (commune a tous les E) ----
FVEC = round(0.5:0.1:20.0, 2);      % 196 frequences, pas 0.1

% ---- Grille k, DEPENDANTE de la bande de frequence ----
%   f <= F_SPLIT (basses freq)          -> k = 0.1..20 (pas 0.1, 200 pts)
%   f >  F_SPLIT (moyennes/hautes freq) -> k = 0.1..60 (pas 0.1, 600 pts)
%   k demarre a 0.1 : k=0 est SINGULIER (termes 1/k dans la matrice A).
F_SPLIT = 2.0;                      % <<< seuil basses | moyennes-hautes (a ajuster)
KLOW  = round(0.1:0.1:20.0, 2);
KHIGH = round(0.1:0.1:60.0, 2);

% ---- TOUS les minima locaux de la courbe Wi_c(k) ----
% -> minima_all_Wi_vs_f_S_0.70.txt : 1 ligne (f, k_c, Wi_c) PAR minimum local.
% MINSEP_K = separation minimale en k entre 2 minima distincts (anti-bruit) :
% baisse-la pour capter + de minima (courbe fine), monte-la pour ne garder que
% les grandes branches.
MINSEP_K = 1.0;

% -------------------- Mode auto-test seul --------------------
% run_S07_matlab('selftest') : verifie le portage et s'arrete (pas de balayage).
if nargin > 0 && (ischar(varargin{1}) || isstring(varargin{1}))
    if strcmpi(varargin{1}, 'selftest')
        selftest_engine(EPS, M_SUB);
        return;
    end
    if strcmpi(varargin{1}, 'crosscheck')
        selftest_engine(EPS, M_SUB);
        fprintf('\nCross-check Wi_c (E=0.006, S=0.7, N=12) vs Python :\n');
        refs = [5.0, 4.073; 10.0, 3.156];   % f, Wi_c_min (Python run_S07_sweep)
        for rr = 1:size(refs,1)
            f  = refs(rr,1);
            kg_ref = 0.5:0.5:40.0;   % grille du run Python (regression du moteur)
            wc = wic_over_k(kg_ref, 2*f^2, 0.006, S, EPS, Nc, M_SUB, WI_LO, 500, 32);
            m  = isfinite(wc) & wc > 0; v = wc(m); kk = kg_ref(m);
            [wm, j] = min(v);
            fprintf('  f=%-4g : MATLAB Wi_c_min=%.4f (k=%.1f)  |  Python %.3f\n', ...
                    f, wm, kk(j), refs(rr,2));
        end
        return;
    end
    if strcmpi(varargin{1}, 'minitest')
        % courbe synthetique : vallees connues a k~5 (locale) et k~15 (globale)
        kv = (0.1:0.1:20)';
        wv = 10 - 4*exp(-((kv-5).^2)/2) - 6*exp(-((kv-15).^2)/2);
        idx = local_minima_idx(wv, kv, 1.0);
        [~, o] = sort(wv(idx)); kc = kv(idx(o));
        fprintf('minitest : %d minima locaux aux k = %s (attendu ~15 puis ~5)\n', ...
                numel(idx), mat2str(round(kc', 1)));
        return;
    end
end

% -------------------- Config par E --------------------
Elist = [0.006, 0.01, 1.0, 10.0, 50.0];
WIHI  = [500,   500,  2000, 12000, 50000];
NSCAN = [32,    32,   36,   40,    44];

% -------------------- Selection des E --------------------
if nargin > 0
    sel = [varargin{:}];
    idx = find(ismember(Elist, sel));
    if isempty(idx), idx = 1:numel(Elist); end
else
    idx = 1:numel(Elist);
end

% -------------------- AUTO-TEST du portage --------------------
selftest_engine(EPS, M_SUB);

% -------------------- Boucle principale --------------------
Ttot = tic;
for ii = idx
    run_one_E(Elist(ii), efmt(Elist(ii)), S, EPS, Nc, M_SUB, ...
              FVEC, KLOW, KHIGH, F_SPLIT, WI_LO, WIHI(ii), NSCAN(ii), MINSEP_K, FORCE);
end
fprintf('\n===== TOUT TERMINE en %.0f s =====\n', toc(Ttot));
end


% =========================================================================
%  DRIVER PAR E
% =========================================================================
function run_one_E(E, tag, S, eps, N, m, fvec, klow, khigh, f_split, wi_lo, wi_hi, nscan, minsep_k, FORCE)
sdir = fullfile('.', ['ucm_E' tag], 'S_0.70');
if ~exist(sdir, 'dir'), mkdir(sdir); end
nf = numel(fvec);
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf(' E=%g  S=%.2f  N=%d  WI_HI=%g  NSCAN=%d\n', E, S, N, wi_hi, nscan);
fprintf(' f=%.2f..%.2f (pas %.2f, %d pts) | k: 0.1..%.0f si f<=%.2f, sinon 0.1..%.0f\n', ...
        fvec(1), fvec(end), fvec(2)-fvec(1), nf, klow(end), f_split, khigh(end));
fprintf('%s\n', repmat('=', 1, 70));

t0 = tic;
parfor ifr = 1:nf
    f    = fvec(ifr);
    fcsv = fullfile(sdir, sprintf('f_%.2f.csv', f));
    if exist(fcsv, 'file') && ~FORCE
        fprintf('  f=%.2f : deja fait (saute)\n', f);
        continue;
    end
    if f <= f_split + 1e-9, kg = klow; else, kg = khigh; end
    wc = wic_over_k(kg, 2.0*f^2, E, S, eps, N, m, wi_lo, wi_hi, nscan);
    write_f_csv(fcsv, kg, wc);
    v = wc(isfinite(wc) & wc > 0);
    if isempty(v)
        fprintf('  f=%.2f : Wi_c_min = NaN   (k<=%.0f)\n', f, kg(end));
    else
        fprintf('  f=%.2f : Wi_c_min = %.4f (k<=%.0f)\n', f, min(v), kg(end));
    end
end

% ----- reconstruit : min GLOBAL (critique) + TOUS les minima locaux -----
crit = fullfile(sdir, 'critical_Wi_vs_f_S_0.70.txt');
mins = fullfile(sdir, 'minima_all_Wi_vs_f_S_0.70.txt');
fid  = fopen(crit, 'w'); fprintf(fid,  'f\tWi_c\tk_c\n');
fid2 = fopen(mins, 'w'); fprintf(fid2, 'f\tk_c\tWi_c\n');   % 1 ligne PAR minimum local
for ifr = 1:nf
    f    = fvec(ifr);
    fcsv = fullfile(sdir, sprintf('f_%.2f.csv', f));
    Wc = NaN; Kc = NaN;
    if exist(fcsv, 'file')
        M  = readmatrix(fcsv);              % colonnes [k, Wi_c]
        kk = M(:,1); ww = M(:,2);
        ok = isfinite(ww) & ww > 0;
        if any(ok)
            kv = kk(ok); wv = ww(ok);
            [Wc, j] = min(wv); Kc = kv(j);          % minimum GLOBAL -> critique
            idxmin = local_minima_idx(wv, kv, minsep_k);   % TOUS les minima locaux
            if ~any(idxmin == j), idxmin = unique([idxmin(:); j]); end
            [~, order] = sort(wv(idxmin));          % du plus critique au moins critique
            for t = order(:)'
                im = idxmin(t);
                fprintf(fid2, '%g\t%g\t%g\n', f, kv(im), wv(im));
            end
        end
    end
    fprintf(fid, '%g\t%g\t%g\n', f, Wc, Kc);
end
fclose(fid); fclose(fid2);
fprintf(' E=%g termine en %.0f s -> %s (+ minima locaux : %s)\n', ...
        E, toc(t0), crit, mins);
end


% =========================================================================
%  Wi_c(k) : rejet parasite + 1er passage neg->pos + fzero bracket
% =========================================================================
function out = wic_over_k(kvec, om, E, S, eps, N, m, wi_lo, wi_hi, nscan)
grid = logspace(log10(wi_lo), log10(wi_hi), nscan);
out  = nan(1, numel(kvec));
for i = 1:numel(kvec)
    k   = kvec(i);
    cst = setup_const(k, om, E, S, eps, N);
    prevW = grid(1);
    prevS = sigma_floq(prevW, cst, k, om, E, S, eps, N, m);
    if prevS >= 0.0
        % sigma(Wi->0) -> negatif pour toute vraie instabilite elastique.
        % sigma(Wi_lo) >= 0 = mode parasite independant de Wi -> rejete (NaN).
        out(i) = NaN;
        continue;
    end
    for jj = 2:nscan
        W = grid(jj);
        s = sigma_floq(W, cst, k, om, E, S, eps, N, m);
        if prevS < 0.0 && s >= 0.0
            try
                out(i) = fzero(@(w) sigma_floq(w, cst, k, om, E, S, eps, N, m), ...
                               [prevW, W]);
            catch
                out(i) = 0.5*(prevW + W);
            end
            break;
        end
        prevW = W; prevS = s;
    end
end
end


% =========================================================================
%  MOTEUR : plus grand exposant de Floquet via monodromie expm renormalisee
% =========================================================================
function sig = sigma_floq(Wi, cst, k, om, E, S, eps, N, m)
% Parties dependantes de Wi (T3..T8), a partir de l'etat de base pre-calcule.
C8 = cst.C8;
T1 = cst.T1; T2 = cst.T2; V1 = cst.V1; V2 = cst.V2; dV1 = cst.dV1; dV2 = cst.dV2;
C9 = Wi; C10 = om*Wi*E;
T3 = -(C9/C8)*T1.*V1 + 2*(C10/C8)*T1.*V2 + 2*(C10/C8)*T2.*V1 + (C9/C8)*T2.*V2;
T4 =  (C9/C8)*T1.*dV1 - 2*(C10/C8)*T1.*dV2 - 2*(C10/C8)*T2.*dV1 - (C9/C8)*T2.*dV2;
T5 = -2*(C10/C8)*T1.*V1 - (C9/C8)*T1.*V2 - (C9/C8)*T2.*V1 - 2*(C10/C8)*T2.*V2;
T6 =  2*(C10/C8)*T1.*dV1 + (C9/C8)*T1.*dV2 + (C9/C8)*T2.*dV1 - 2*(C10/C8)*T2.*dV2;
T7 = -C9*(T1.*V1 + T2.*V2);
T8 =  C9*(T1.*dV1 + T2.*dV2);

D1 = cst.D1; D2 = cst.D2; D4 = cst.D4; I = cst.I; Z = cst.Z;
G = cst.G; invBB = cst.invBB; r = cst.r; ki = cst.ki;

T  = 2*pi/om; dt = T/m;
Phi = eye(cst.nred);
logscale = 0.0;
for j = 1:m
    t  = (j - 0.5)*dt;
    Vi = V1*cos(om*t) + V2*sin(om*t);
    Ti = T1*cos(om*t) + T2*sin(om*t);
    Tii = (eps*T3 + T4)*cos(2*om*t) + (eps*T5 + T6)*sin(2*om*t) + (eps*T7 + T8);
    Vc = diag(Vi); dVs = diag(D1*Vi); Tc = diag(Ti); dTs = diag(D1*Ti);
    Tcc = diag(Tii); dTss = diag(D1*Tii);

    A = complex(zeros(8*N));
    A(3:N-2,:) = [(1-S)*(D4(3:N-2,:) - 2*k^2*D2(3:N-2,:) + k^4*I(3:N-2,:)), Z(3:N-2,:), -k^2*D1(3:N-2,:), Z(3:N-2,:), -1i*k^3*I(3:N-2,:) - 1i*k*D2(3:N-2,:), eps*k^2*I(3:N-2,:), Z(3:N-2,:), k^2*D1(3:N-2,:)];
    A(N+2:2*N-1,:) = [Z(2:N-1,1:N), (1-S)*(D2(2:N-1,:) - k^2*I(2:N-1,:)), Z(2:N-1,:), D1(2:N-1,:) + 2*eps*I(2:N-1,:), Z(2:N-1,:), Z(2:N-1,:), 1i*k*I(2:N-1,:), Z(2:N-1,:)];
    A(2*N+1:3*N,:) = [2*S*D1, Z, -I, Z, Z, Z, Z, Z];
    A(3*N+1:4*N,:) = [Wi*Tc*D1 - Wi*dTs + eps*Wi*Tc, S*D1, Wi*dVs - eps*Wi*Vc, -I, Z, Z, Z, Z];
    A(4*N+1:5*N,:) = [1i*(1/k)*S*D2 + 1i*k*S*I, Z, Z, Z, -I, Z, Z, Z];
    A(5*N+1:6*N,:) = [-Wi*dTss + 2*eps*Wi*Tcc, 2*Wi*Tc*D1 - 2*eps*Wi*Tc, Z, 2*Wi*dVs - 2*eps*Wi*Vc, Z, -I, Z, Z];
    A(6*N+1:7*N,:) = [1i*(1/k)*Wi*(Tc*D2 + eps*Tc*D1 - eps^2*Tc), 1i*k*S*I, Z, Z, Wi*dVs - eps*Wi*Vc, Z, -I, Z];
    A(7*N+1:8*N,:) = [-2*S*D1 - 2*eps*S*I, Z, Z, Z, Z, Z, Z, -I];

    AA  = A(ki,ki) + A(ki,r)*G;
    qqq = AA*invBB;
    Phi = expm(qqq*dt)*Phi;
    nrm = max(abs(Phi(:)));
    if isfinite(nrm) && nrm > 0
        Phi = Phi/nrm;
        logscale = logscale + log(nrm);
    end
end
if ~all(isfinite(Phi(:)))
    sig = 1e3; return;      % croissance enorme -> tres instable
end
mu  = eig(Phi);
sig = (max(log(abs(mu))) + logscale)/T;
end


% =========================================================================
%  Pre-calcul des parties INDEPENDANTES de Wi et t (une fois par k,om,E,S)
% =========================================================================
function cst = setup_const(k, om, E, S, eps, N)
[x, DM] = chebdif(N, 4);
s = 2;
D1 = s*DM(:,:,1); D2 = s^2*DM(:,:,2); D4 = s^4*DM(:,:,4);
I = eye(N); Z = zeros(N);

% ----- etat de base (verbatim rhs_B) -----
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

% ----- matrice B (constante) -----
B = zeros(8*N);
B(3:N-2,:) = [D2(3:N-2,1:N) - k^2*I(3:N-2,1:N), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:), Z(3:N-2,:)];
B(N+2:2*N-1,:) = [Z(2:N-1,1:N), I(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N), Z(2:N-1,1:N)];
for blk = 3:8
    B((blk-1)*N+1:blk*N, (blk-1)*N+1:blk*N) = E*I;
end

% ----- conditions aux limites / reduction -----
c1 = [I(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:)];
c2 = [D1(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:)];
c3 = [I(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:)];
c4 = [D1(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:)];
c5 = [Z(1,:), I(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:), Z(1,:)];
c6 = [Z(N,:), I(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:), Z(N,:)];
Cmat = [c1; c2; c3; c4; c5; c6];
r  = [1, 2, N-1, N, N+1, 2*N];
ki = [3:N-2, N+2:2*N-1, 2*N+1:8*N];
G  = -pinv(Cmat(:,r))*Cmat(:,ki);
BB = B(ki,ki) + B(ki,r)*G;

cst = struct('D1',D1,'D2',D2,'D4',D4,'I',I,'Z',Z, ...
             'V1',V1,'V2',V2,'dV1',dV1,'dV2',dV2,'T1',T1,'T2',T2,'C8',C8, ...
             'r',r,'ki',ki,'G',G,'invBB',inv(BB),'nred',8*N-6);
end


% =========================================================================
%  Matrices de differentiation Chebyshev-Gauss-Lobatto (Weideman-Reddy)
% =========================================================================
function [x, DM] = chebdif(N, M)
I = eye(N); L = logical(I);
n1 = floor(N/2); n2 = ceil(N/2);
k = (0:N-1)'; th = k*pi/(N-1);
x = sin(pi*(N-1:-2:1-N)'/(2*(N-1)));
T = repmat(th/2, 1, N);
DX = 2*sin(T'+T).*sin(T'-T);
DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))]; DX(L) = ones(N,1);
C = toeplitz((-1).^k);
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2; C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;
Zc = 1./DX; Zc(L) = zeros(N,1);
D = eye(N);
for ell = 1:M
    D = ell*Zc.*(C.*repmat(diag(D),1,N) - D);
    D(L) = -sum(D');
    DM(:,:,ell) = D;
end
end


% =========================================================================
%  Auto-test : sigma du portage vs references du moteur Python valide
% =========================================================================
function selftest_engine(eps, m)
fprintf('Auto-test du moteur (sigma vs references Python valide)...\n');
cases = {
%   E,     S,    om,        k,    Wi,     N,   sigma_ref
    0.006, 0.7,  2*5.0^2,   6.0,  3.6,    8,   -17.4466691934
    1.0,   1.0,  2*2.0^2,   14.0, 4.7,    12,   1.7123623417
    10.0,  0.7,  2*2.0^2,   22.0, 1000.0, 12,  -0.0399338047
};
okall = true;
for c = 1:size(cases,1)
    E = cases{c,1}; S = cases{c,2}; om = cases{c,3};
    k = cases{c,4}; Wi = cases{c,5}; N = cases{c,6}; ref = cases{c,7};
    cst = setup_const(k, om, E, S, eps, N);
    s   = sigma_floq(Wi, cst, k, om, E, S, eps, N, m);
    err = abs(s - ref);
    tol = 1e-3 + 1e-4*abs(ref);
    stat = 'OK'; if err > tol, stat = 'ECHEC'; okall = false; end
    fprintf('  E=%-5g S=%.2f N=%-2d Wi=%-7g : sigma=%.8f (ref %.8f, err %.2e) [%s]\n', ...
            E, S, N, Wi, s, ref, err, stat);
end
if ~okall
    error('run_S07_matlab:selftest', ...
          'AUTO-TEST ECHOUE : le portage ne reproduit pas le moteur valide. Abandon.');
end
fprintf('Auto-test OK.\n');
end


% =========================================================================
%  Utilitaires
% =========================================================================
function tag = efmt(E)
tag = strrep(strrep(sprintf('%g', E), '.', 'p'), '-', 'm');
end

function write_f_csv(fname, kgrid, wc)
fid = fopen(fname, 'w');
fprintf(fid, 'k, Wi_c\n');
for ik = 1:numel(kgrid)
    fprintf(fid, '%g, %g\n', kgrid(ik), wc(ik));
end
fclose(fid);
end

function idx = local_minima_idx(wv, kv, sep)
% Indices des minima LOCAUX de la courbe wv (=Wi_c) en fonction de kv (=k),
% separes d'au moins 'sep' en k (islocalmin garde le plus profond par fenetre).
wv = wv(:); kv = kv(:);
if numel(wv) < 3
    [~, idx] = min(wv);
    return;
end
try
    tf = islocalmin(wv, 'SamplePoints', kv, 'MinSeparation', sep);
catch
    % repli si islocalmin indisponible : minima locaux stricts (interieurs)
    tf = [false; (wv(2:end-1) < wv(1:end-2)) & (wv(2:end-1) < wv(3:end)); false];
end
idx = find(tf);
if isempty(idx)
    [~, idx] = min(wv);     % au minimum, le min global
end
end
