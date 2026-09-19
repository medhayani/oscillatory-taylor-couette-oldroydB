function sigma = sigma_B_floquet(Wi, k, om, E, S, eps)
% =========================================================================
% Calcul du plus grand exposant de Floquet (partie reelle) pour le systeme B
% (formulation visqueuse corrigee avec terme de courbure elastique)
%
% INPUTS :
%   Wi    : Weissenberg de controle
%   k     : nombre d'onde axial q
%   om    : frequence visqueuse sigma = omega*d^2/nu = 2*gamma^2
%   E     : nombre elastique = lambda*nu/d^2
%   S     : rapport polymere/total = eta_p/eta
%   eps   : rapport de gap d/R_1
%
% OUTPUT :
%   sigma : max_j Re(sigma_j) ou sigma_j = log(mu_j)/T sont les exposants
%           de Floquet, mu_j etant les valeurs propres de la matrice de
%           monodromie sur la periode T = 2*pi/om
% =========================================================================

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
% RHS de l'equation de la matrice fondamentale d/dt phi = qqq * phi
% avec qqq = AA / BB (apres elimination des CL et reduction matricielle)

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
