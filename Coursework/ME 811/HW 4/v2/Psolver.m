function [pp Resp] = Psolver(ApO,ApE,ApW,ApN,ApS,Sp)
    Rtol = 1E-20;
    itrmax = 20;
    [M N] = size(Sp);
%     %Reshape link coefficients
%     ApO = reshape(ApO,N*M,1);
%     ApE = reshape(ApE,N*M,1);
%     ApW = reshape(ApW,N*M,1);
%     ApN = reshape(ApN,N*M,1);
%     ApS = reshape(ApS,N*M,1);
%     
%     ApN = circshift(ApN,-1);
%     ApS = circshift(ApS,1);%need to check to make sure AN,AS are being shifted properly
%     ApE = circshift(ApE,M);
%     ApW = circshift(ApW,-M);
%     
%     Ap = spdiags([ApW ApN ApO ApS ApE],[-M -1 0 1 M],N*M,N*M);
%     
%     %Solve for u hat
%     [pp R] = ADI2d(Ap,Sp,'-TDMA',Rtol,itrmax);
%     Rp = norm(R);
    ApO = flipud(ApO);
    ApE =  flipud(ApE);
    ApW =  flipud(ApW);
    ApN =  flipud(ApN);
    ApS =  flipud(ApS);
    Sp = flipud(Sp);
%     [pp R] = ADIp(ApO,ApE,ApW,ApN,ApS,Sp,Rtol,itrmax);
    [pp,R] = ADI(ApS,ApW,ApO,ApE,ApN,Sp,zeros(M,N),Rtol,itrmax);
    pp = flipud(pp);
    Resp = R(end);
end