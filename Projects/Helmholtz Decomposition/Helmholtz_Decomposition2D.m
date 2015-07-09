function [potential,solenoidal] = Helmholtz_Decomposition2D(x,y,U,V)

    N = size(U);
    
    %Make sure matrices are in correct ascending order
    [~,Ix] = sort(x(:,1));
    [~,Iy] = sort(y(1,:));
    x = x(Ix,:);
    y = y(:,Iy);
    U = U(Ix,Iy);
    V = V(Ix,Iy);
    
    %Create derivative matrices for divergence calculation
    order = 4;
    dx = mean(diff(x(:,1)));
    dy = mean(diff(y(1,:))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    partial_x_D1 = mNumericalDerivative(1,order,dx,N(1));
    partial_y_D1 = mNumericalDerivative(1,order,dy,N(2));
    
    divergence = partial_x_D1*U + (partial_y_D1*V')'; %for now we are just going to assume that U and V are 2-D
    
    %Create derivative matrix for poisson equation, A*Phi = F
    L = ceil(order/2);
    NM = prod(N(1:2));
    partial_x_D2_coefs = TSE(L,L,dx,2);
    partial_y_D2_coefs = TSE(L,L,dy,2);
    Ax = spdiags(repmat(partial_x_D2_coefs(:)',NM,1),-L:L, NM,NM);
    Ay = spdiags(repmat(partial_y_D2_coefs(:)',NM,1),(-L:L)*N(1), NM,NM); %offset on second group diagonals is due to reshape of N X M matrix - note that we still need to fix the boundaries
    F = divergence(:);
    %Now we need to fix the boundary nodes...
    %left boundary
    partial_x_D1_left = TSE(0,2*L,dx,1);
    
    
end