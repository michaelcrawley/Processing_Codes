function [potential,solenoidal] = Helmholtz_Decomposition2DCyl_v3(z,r,Uz,Ur)
    %This function assumes that the lower boundary is the jet centerline,
    %and that the potential flow component is zero at the inlet and outlet
    %boundaries.

    [M,N] = size(Uz);
    %First DIM should be y, second DIM x
    
    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(z(1,:));
    [~,Ir] = sort(r(:,1));
    z = z(:,Iz);
    r = r(Ir,:);
    Uz = Uz(Ir,Iz);
    Ur = Ur(Ir,Iz);
        
    %Create derivative matrices for divergence calculation
    order = 4;
    dz = mean(diff(z(1,:)));
    dr = mean(diff(r(:,1))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    partial_z_D1 = mNumericalDerivative(1,order,dz,N);
    partial_r_D1 = mNumericalDerivative(1,order,dr,M);
    
    divergence = partial_r_D1*Ur + (partial_z_D1*Uz')' + Ur./r; %for now we are just going to assume that U and V are 2-D
    %Get integrate derivatives at upper and lower boundaries for boundary
    %conditions
    int_top = cumsum(Uz(end,:))*dz;
    
    %Create derivative matrix for poisson equation, A*Phi = F
    L = 1; %for now, we're just going to use a second order method for the poisson solver 
    zcoefs = [1 -2 1]/dz/dz;
    rcoefs1 = [1 0 -1]/2/dr; %we still need to add in the 1/r dependence on these coefficients
    rcoefs = [1 -2 1]/dr/dr;
    Az = spdiags(repmat(zcoefs,N*M,1),(-L:L)*M, N*M,N*M);
    Ar = spdiags(repmat(rcoefs,N*M,1),(-L:L), N*M,N*M) + spdiags(repmat(rcoefs1,N*M,1)./repmat(r(:),1,3),(-L:L), N*M,N*M)'; %offset on second group diagonals is due to reshape of N X M matrix - note that we still need to fix the boundaries
    F = divergence(:);
    A = Az + Ar;
    
    %Now we need to fix the boundary nodes...
    %Bottom - (1,:)
%     A(1:M:1+(N-1)*M,1:M:1+(N-1)*M) = spdiags(repmat(zcoefs,N,1),(-L:L)*M, N,N) + spdiags(repmat([2 -2 0]/dr/dr,N,1),(-L:L), N,N);
    for n = 1+M:M:1+(N-1)*M
        A(n,n-1) = 0;
        A(n,n+1) = 2/dr/dr;
    end    
    %Top - (M,:)
    A(M:M:M+(N-1)*M,:) = 0;
    A(M:M:M+(N-1)*M,M:M:M+(N-1)*M) = sparse(eye(N));
    F(M:M:M+(N-1)*M) = int_top;
    %Inflow - (:,1)
    A(1,1:M+1) = sparse([-1/dr-1/dz,1/dr,zeros(1,M-2),1/dz]);
    for n = 2:M-1
        A(n,n-1:M+n) = [-1/2/dr,-1/dz,1/2/dr,zeros(1,M-2),1/dz];
    end
    A(M,M-1:2*M) = [-1/dr,-1/dr-1/dz,zeros(1,M-1),1/dz];
    F(1:M) = 0;
    %Outflow - (:,N)
    br_corner = 1+(N-1)*M;
    tr_corner = N*M;
    A(br_corner,br_corner-M:br_corner+1) = [-1/dz,zeros(1,M-1),1/dz-1/dr,1/dr];
    for n = br_corner+1:tr_corner-1
        A(n,n-M:n+1) = [-1/dz,zeros(1,M-2),-1/2/dr,1/dz,1/2/dr];
    end
    A(tr_corner,tr_corner-M:tr_corner) = [-1/dz,zeros(1,M-2),-1/dr,1/dr+1/dz];
    F(end-M+1:end) = 0;
    
    %
    phi = reshape(A\F,M,N);
    potential.Ur = partial_r_D1*phi;
    potential.Uz = (partial_z_D1*phi')';
    solenoidal.Ur = Ur - potential.Ur;
    solenoidal.Uz = Uz - potential.Uz;
end