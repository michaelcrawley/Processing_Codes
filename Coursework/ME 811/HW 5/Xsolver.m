function [uh Resx] = Xsolver(alpha,AxO,AxE,AxW,AxN,AxS,Sx,u)
    Rtol = 1E-20;
    itrmax = 2;
    [M N] = size(Sx);
    %Reshape link coefficients
    AxO = reshape(AxO,N*M,1);
    AxE = reshape(AxE,N*M,1);
    AxW = reshape(AxW,N*M,1);
    AxN = reshape(AxN,N*M,1);
    AxS = reshape(AxS,N*M,1);
    
    AxN = circshift(AxN,-1);
    AxS = circshift(AxS,1);%need to check to make sure AN,AS are being shifted properly
    AxE = circshift(AxE,M);
    AxW = circshift(AxW,-M);
    
    Ax = spdiags([AxW AxN AxO AxS AxE],[-M -1 0 1 M],N*M,N*M);
    
    %Modify source for correction equation
    Rx = reshape(Sx,M*N,1)-Ax*reshape(u(:,2:end),M*N,1);
    Sx = reshape(Rx,M,N);
    Resx = norm(Rx);
    Ax = Ax+diag(AxO*alpha);
    
    %Solve for u hat
    [uh R] = ADIp(Ax,Sx,'-TDMA',Rtol,itrmax);

    %append boundary condition
    uh = [zeros(M,1) uh];
end

