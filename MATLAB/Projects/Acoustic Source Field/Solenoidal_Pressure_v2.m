function [p_s] = Solenoidal_Pressure_v2(z,r,rho,s_Uz,s_Ur)   
%This function assumes that the lower boundary is the jet centerline,
%and the solenoidal pressure is zero at the rest of the boundaries.

    [M,N] = size(s_Uz);
    %First DIM should be y, second DIM x
    
    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(z(1,:));
    [~,Ir] = sort(r(:,1));
    z = z(:,Iz);
    r = r(Ir,:);
    s_Uz = s_Uz(Ir,Iz);
    s_Ur = s_Ur(Ir,Iz);
    dz = mean(diff(z(1,:)));
    dr = mean(diff(r(:,1))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    
    %Compute incompressible source field
    S = ComputeAASource(r,z,rho,s_Ur,s_Uz);
    
    %Create derivative matrix for poisson equation, A*Phi = F
    L = 1; %for now, we're just going to use a second order method for the poisson solver 
    zcoefs = [1 -2 1]/dz/dz;
    rcoefs1 = [1 0 -1]/2/dr; %we still need to add in the 1/r dependence on these coefficients
    rcoefs = [1 -2 1]/dr/dr;
    Az = spdiags(repmat(zcoefs,N*M,1),(-L:L)*M, N*M,N*M);
    Ar = spdiags(repmat(rcoefs,N*M,1),(-L:L), N*M,N*M) + spdiags(repmat(rcoefs1,N*M,1)./repmat(r(:),1,3),(-L:L), N*M,N*M)'; %offset on second group diagonals is due to reshape of N X M matrix - note that we still need to fix the boundaries
    F = -S(:);
    A = Az + Ar;
    
    %Now we need to fix the boundary nodes...
    %Bottom - (1,:)
    for n = 1+M:M:1+(N-1)*M
        A(n,n-1) = 0;
        A(n,n+1) = 2/dr/dr;
    end    
    %Top - (M,:)
    A(M:M:M+(N-1)*M,:) = 0;
    A(M:M:M+(N-1)*M,M:M:M+(N-1)*M) = sparse(eye(N));
    F(M:M:M+(N-1)*M) = 0;
    %Inflow - (:,1)
%     A(1,1:M+1) = [-2/dz/dz-2/dr/dr,2/dr/dr,zeros(1,M-2),2/dz/dz];
%     for n = 2:M-1
%         A(n,n+M) = 2/dz/dz;
%     end
    A(1:M,:) = 0;
    A(1:M,1:M) = eye(M);
    F(1:M) = 0;
    %Outflow - (:,N)
    for n = (N-1)*M+2:N*M-1
        A(n,n-M) = 2/dz/dz;
    end
    corner = 1+(N-1)*M;
    A(corner,corner-M:corner+1) = [2/dz/dz,zeros(1,M-1),-2/dz/dz-2/dr/dr,2/dr/dr];
    
    %
    p_s = reshape(A\F,M,N);    
end