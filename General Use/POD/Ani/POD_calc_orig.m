function [Lambda, Phi] = POD_calc_orig(u, inner_product, W)
%W is the weighting matrix for the inner product, but since it is assumed
%to be diagonal, it is also assumed to be passed to this function as a
%vector. Also, the same matrix is passed on to the inner_product routine
%too. See Rowley's thesis eqn 4.9 for an explanation of how the same
%weighting matrix applies to the POD problem.

[N,M] = size(u);

switch nargin
    case 2
        W_half = eye(N);
        W_half_inv = eye(N);
    case 3
        sqrtW = sqrt(W);
        W_half = diag(sqrtW);
        W_half_inv = diag(1./sqrtW);
    otherwise
        error('Incorrect number of arguments')
end

%Correlation matrix is simply the outer-product divided by the number of
%snapshots
R = (u*u.')/M;
RW = W_half*R*W_half; %Weighting matrix applied per Glauser's formulation

[U,S,V] = svd(RW);
Lambda = diag(S);
Phi = W_half_inv*V; %Revert back to physical eigenvectors

%Normalize Phi by their 2-norm
Phi_2_norm = sqrt(diag(inner_product(Phi,Phi,W)));
Phi = Phi*diag(1./Phi_2_norm);
