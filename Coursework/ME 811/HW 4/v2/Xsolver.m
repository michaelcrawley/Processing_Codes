function [uh Resx] = Xsolver(alpha,AxO,AxE,AxW,AxN,AxS,Sx,u)
    Rtol = 1E-20;
    itrmax = 2;
    [M N] = size(Sx);
%     %Reshape link coefficients
%     AxO = reshape(AxO,N*M,1);
%     AxE = reshape(AxE,N*M,1);
%     AxW = reshape(AxW,N*M,1);
%     AxN = reshape(AxN,N*M,1);
%     AxS = reshape(AxS,N*M,1);
%     
%     AxN = circshift(AxN,-1);
%     AxS = circshift(AxS,1);%need to check to make sure AN,AS are being shifted properly
%     AxE = circshift(AxE,M);
%     AxW = circshift(AxW,-M);
%     
%     Ax = spdiags([AxW AxN AxO AxS AxE],[-M -1 0 1 M],N*M,N*M);
%     
%     %Modify source for correction equation
%     Rx = reshape(Sx,M*N,1)-Ax*reshape(u(:,2:end),M*N,1);
%     Sx = reshape(Rx,M,N);
%     Resx = norm(Rx);
%     Ax = Ax+diag(AxO*alpha);
%     
%     %Solve for u hat
%     [uh R] = ADI2d(Ax,Sx,'-TDMA',Rtol,itrmax);
    AxO = flipud([ones(M,1) AxO ones(M,1)]);
    AxE = flipud([zeros(M,1) AxE zeros(M,1)]);
    AxW = flipud([zeros(M,1) AxW zeros(M,1)]);
    AxS = flipud([zeros(M,1) AxS zeros(M,1)]);
    AxN = flipud([zeros(M,1) AxN zeros(M,1)]);
    Sx = flipud([zeros(M,1) -Sx zeros(M,1)]);
    u = flipud([u zeros(M,1)]);
   
    
    [Resx,Ri] = calcRes(AxO,AxE,AxW,AxN,AxS,Sx,u);
%     [uh,~] = ADIp((1+alpha)*AxO,AxE,AxW,AxN,AxS,Ri,Rtol,itrmax);
   [uh,~] = ADI(AxS,AxW,(1+alpha)*AxO,AxE,AxN,Ri,zeros(M,N+2),Rtol,itrmax);
    
    %append boundary condition
    uh = flipud(uh(:,1:end-1));
end

