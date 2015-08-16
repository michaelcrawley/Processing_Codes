function [potential,solenoidal] = Helmholtz_Decomposition2DCart(x,y,U,V)

    [M,N] = size(U);
    %First DIM should be y, second DIM x
    
    %Make sure matrices are in correct ascending order
    [~,Ix] = sort(x(1,:));
    [~,Iy] = sort(y(:,1));
    x = x(:,Ix);
    y = y(Iy,:);
    U = U(Iy,Ix);
    V = V(Iy,Ix);
        
    %Create derivative matrices for divergence calculation
    order = 4;
    dx = mean(diff(x(1,:)));
    dy = mean(diff(y(:,1))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    partial_x_D1 = mNumericalDerivative(1,order,dx,N);
    partial_y_D1 = mNumericalDerivative(1,order,dy,M);
    
    divergence = partial_y_D1*V + (partial_x_D1*U')'; %for now we are just going to assume that U and V are 2-D
    %Get integrate derivatives at upper and lower boundaries for boundary
    %conditions
    int_top = cumsum(U(1,:))*dx;
    int_bottom = cumsum(U(end,:))*dx;
    
    %Create derivative matrix for poisson equation, A*Phi = F
    L = 1; %for now, we're just going to use a second order method for the poisson solver 
    xcoefs = [1 -2 1]/dx/dx;
    ycoefs = [1 -2 1]/dy/dy;
    Ax = spdiags(repmat(xcoefs,N*M,1),(-L:L)*M, N*M,N*M);
    Ay = spdiags(repmat(ycoefs,N*M,1),(-L:L), N*M,N*M); %offset on second group diagonals is due to reshape of N X M matrix - note that we still need to fix the boundaries
    F = divergence(:);
    A = Ax + Ay;
    
    %Now we need to fix the boundary nodes...
    %Top - (1,:)
    A(1:M:1+(N-1)*M,1:M:1+(N-1)*M) = sparse(eye(N));
    F(1:M:1+(N-1)*M) = int_top;
    %Bottom - (M,:)
    A(M:M:M+(N-1)*M,M:M:M+(N-1)*M) = sparse(eye(N));
    F(M:M:M+(N-1)*M) = int_bottom;
    %Inflow - (:,1)
    A((2:M-1),M+(2:M-1)) = 2*A((2:M-1),M+(2:M-1));
    F(2:M-1) = F(2:M-1) + 2*dx*U(2:M-1,1);
    %Outflow - (:,N)
    A((N-1)*M + (2:M-1),(N-2)*M +(2:M-1)) = 2*A((N-1)*M + (2:M-1),(N-2)*M +(2:M-1));
    F((N-1)*M+2:(N-1)*M+M-1) = F((N-1)*M+2:(N-1)*M+M-1) + 2*dx*U(2:M-1,N);
    
    %
    phi = reshape(A\F,M,N);
    potential.v = partial_y_D1*phi;
    potential.u = (partial_x_D1*phi')';
    solenoidal.v = V - potential.v;
    solenoidal.u = U - potential.u;
end