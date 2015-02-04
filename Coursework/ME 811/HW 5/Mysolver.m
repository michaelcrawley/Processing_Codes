function [uh Resx] = Mysolver(alpha,AO,AE,AW,AN,AS,S,u)
    Rtol = 1E-20;
    itrmax = 2;
    [M N] = size(S);
    %Reshape link coefficients
    AOt = reshape(AO,N*M,1);
    AEt = reshape(AE,N*M,1);
    AWt = reshape(AW,N*M,1);
    ANt = reshape(AN,N*M,1);
    ASt = reshape(AS,N*M,1);
    
    ANt = circshift(ANt,-1);
    ASt = circshift(ASt,1);%need to check to make sure AN,AS are being shifted properly
    AEt = circshift(AEt,M);
    AWt = circshift(AWt,-M);
    
    Ax = spdiags([AWt ANt AOt ASt AEt],[-M -1 0 1 M],N*M,N*M);
    
    %Modify source for correction equation
    Rx = reshape(S,M*N,1)-Ax*reshape(u(1:end-1,:),M*N,1);
    Sx = reshape(Rx,M,N);
    Resx = norm(Rx);
    
    %Solve for u hat
%     [uh R] = ADIp(Ax,Sx,'-TDMA',Rtol,itrmax);
    uh = ADIc( AN , AW , (1+alpha)*AO , AE , AS , Rx , zeros(M,N) , -inf , itrmax );

    %append boundary condition
    uh = [zeros(M,1) uh];
end