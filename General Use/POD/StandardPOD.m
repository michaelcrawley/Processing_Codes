function [phi,lambda,ak] = StandardPOD(U,wk)
%DIM1 of U is position
%DIM2 of U is time
%wk is a weight vector, which can be used for non-linearly spaced positions
%or times

    N = size(U);
    if N(2) < N(1), warning('Either U matrix should be rotated, or SnapShotPOD method should be used'); end
    if ~exist('wk','var')||isempty(wk), wk = ones(N); end
    M = size(wk);
    if any(M == 1), wk = repmat(wk,N./M); end
    U = U.*sqrt(wk); %energy scaling, sqrt as it is applied before R matrix
    
    R = (U*U')/N(2); %cross-correlation matrix over time - note that conjugate transpose is intentional
    [phi,eigs] = svd(R);    lambda = diag(eigs);
    ak = (U.'*phi).';
end