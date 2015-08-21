function [p_s] = Solenoidal_Pressure(z,r,rho,Uz,Ur)   
%This function assumes that the lower boundary is the jet centerline,
%and the solenoidal pressure is zero at the rest of the boundaries.

    [M,N] = size(Uz);
    %First DIM should be y, second DIM x
    
    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(z(1,:));
    [~,Ir] = sort(r(:,1));
    z = z(:,Iz);
    r = r(Ir,:);
    Uz = Uz(Ir,Iz);
    Ur = Ur(Ir,Iz);
    dz = mean(diff(z(1,:)));
    dr = mean(diff(r(:,1))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    
    %Compute incompressible source field
    S = ComputeAASource(r,z,rho,Ur,Uz);
    
    %Create derivative matrix for poisson equation, A*Phi = F
    L = 1; %for now, we're just going to use a second order method for the poisson solver 
    zcoefs = [1 -2 1]/dz/dz;
    rcoefs1 = [-1 0 1]/2/dr; %we still need to add in the 1/r dependence on these coefficients
    rcoefs = [1 -2 1]/dr/dr;
    Az = spdiags(repmat(zcoefs,N*M,1),(-L:L)*M, N*M,N*M);
    Ar = spdiags(repmat(rcoefs,N*M,1),(-L:L), N*M,N*M) + spdiags(repmat(rcoefs1,N*M,1)./repmat(r(:),1,3),(-L:L), N*M,N*M); %offset on second group diagonals is due to reshape of N X M matrix - note that we still need to fix the boundaries
    F = -S(:);
    A = Az + Ar;
    
    %Now we need to fix the boundary nodes...
    %Bottom - (1,:)
    A(1:M:1+(N-1)*M,1:M:1+(N-1)*M) = spdiags(repmat(zcoefs,N,1),(-L:L)*M, N,N) + spdiags(repmat([2 -2 0]/dr/dr,N,1),(-L:L), N,N);
    %Top - (M,:)
    A(M:M:M+(N-1)*M,M:M:M+(N-1)*M) = sparse(eye(N));
    F(M:M:M+(N-1)*M) = 0;
    %Inflow - (:,1)
    A((2:M-1),(2:M-1)) = sparse(eye(M-2));
    F(2:M-1) = 0;
    %Outflow - (:,N)
    A(end-(2:M-1),end-(2:M-1)) = sparse(eye(M-2));
    F((N-1)*M+2:(N-1)*M+M-1) = 0;
    
    %
    p_s = reshape(A\F,M,N);    
end