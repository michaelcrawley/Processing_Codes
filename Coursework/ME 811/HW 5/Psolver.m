function [pp Rp] = Psolver(ApO,ApE,ApW,ApN,ApS,Sp)
    Rtol = 1E-20;
    itrmax = 20;
    [M N] = size(Sp);
    %Reshape link coefficients
    ApO = reshape(ApO,N*M,1);
    ApE = reshape(ApE,N*M,1);
    ApW = reshape(ApW,N*M,1);
    ApN = reshape(ApN,N*M,1);
    ApS = reshape(ApS,N*M,1);
    
    ApN = circshift(ApN,-1);
    ApS = circshift(ApS,1);%need to check to make sure AN,AS are being shifted properly
    ApE = circshift(ApE,M);
    ApW = circshift(ApW,-M);
    
    Ap = spdiags([ApW ApN ApO ApS ApE],[-M -1 0 1 M],N*M,N*M);
    
    %Solve for u hat
    [pp R] = ADIp(Ap,Sp,'-TDMA',Rtol,itrmax);
    Rp = norm(R);
end