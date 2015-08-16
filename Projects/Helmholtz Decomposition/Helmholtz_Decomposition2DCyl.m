function [potential,solenoidal] = Helmholtz_Decomposition2DCyl(z,r,Uz,Ur)

    [M,N] = size(Uz);
    %First DIM should be y, second DIM x
    
    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(z(1,:));
    [~,Ir] = sort(r(:,1));
    z = z(:,Iz);
    r = r(Ir,:);
    Uz = Uz(Ir,Iz);
    Ur = Ur(Ir,Iz);
    Ur(r<0) = -1*(Ur(r<0)); %there really shouldn't be negative radial points, so we need to flip the signs where we do
        
    %Create derivative matrices for divergence calculation
    order = 4;
    dz = mean(diff(z(1,:)));
    dr = mean(diff(r(:,1))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    partial_z_D1 = mNumericalDerivative(1,order,dz,N);
    partial_r_D1 = mNumericalDerivative(1,order,dr,M);
    r = abs(r); %really, r shouldn't be negative, its just cuz we don't have theta that we are using this convention
    
    divergence = partial_r_D1*Ur + (partial_z_D1*Uz')' + Ur./r; %for now we are just going to assume that U and V are 2-D
    %Get integrate derivatives at upper and lower boundaries for boundary
    %conditions
    int_bottom = cumsum(Uz(1,:))*dz;
    int_top = cumsum(Uz(end,:))*dz;
    
    %Create derivative matrix for poisson equation, A*Phi = F
    L = 1; %for now, we're just going to use a second order method for the poisson solver 
    zcoefs = [1 -2 1]/dz/dz;
    rcoefs1 = [-1 0 1]/2/dr; %we still need to add in the 1/r dependence on these coefficients
    rcoefs = [1 -2 1]/dr/dr;
    Az = spdiags(repmat(zcoefs,N*M,1),(-L:L)*M, N*M,N*M);
    Ar = spdiags(repmat(rcoefs,N*M,1),(-L:L), N*M,N*M) + spdiags(repmat(rcoefs1,N*M,1)./repmat(r(:),1,3),(-L:L), N*M,N*M); %offset on second group diagonals is due to reshape of N X M matrix - note that we still need to fix the boundaries
    F = divergence(:);
    A = Az + Ar;
    
    %Now we need to fix the boundary nodes...
    %Bottom - (1,:)
    A(1:M:1+(N-1)*M,1:M:1+(N-1)*M) = sparse(eye(N));
    F(1:M:1+(N-1)*M) = int_bottom;
    %Top - (M,:)
    A(M:M:M+(N-1)*M,M:M:M+(N-1)*M) = sparse(eye(N));
    F(M:M:M+(N-1)*M) = int_top;
    %Inflow - (:,1)
    A((2:M-1),M+(2:M-1)) = 2*A((2:M-1),M+(2:M-1));
    F(2:M-1) = F(2:M-1) + 2*Uz(2:M-1,1)/dz;
    %Outflow - (:,N)
    A((N-1)*M + (2:M-1),(N-2)*M +(2:M-1)) = 2*A((N-1)*M + (2:M-1),(N-2)*M +(2:M-1));
    F((N-1)*M+2:(N-1)*M+M-1) = F((N-1)*M+2:(N-1)*M+M-1) - 2*Uz(2:M-1,N)/dz;
    
    %
    phi = reshape(A\F,M,N);
    potential.v = partial_r_D1*phi;
    potential.u = (partial_z_D1*phi')';
    solenoidal.v = Ur - potential.v;
    solenoidal.u = Uz - potential.u;
end