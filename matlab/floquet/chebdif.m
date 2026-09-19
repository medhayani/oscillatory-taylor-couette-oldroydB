function [x, DM] = chebdif(N, M)
% =========================================================================
% Matrices de differentiation Chebyshev-Gauss-Lobatto (Weideman-Reddy)
%
% INPUTS :
%   N : nombre de points de collocation
%   M : ordre maximum des derivees a calculer
%
% OUTPUTS :
%   x  : vecteur des points de Chebyshev sur [-1, 1]
%   DM : array 3D de matrices de differentiation, DM(:,:,ell) = matrice
%        d'ordre ell pour ell = 1, ..., M
% =========================================================================

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
