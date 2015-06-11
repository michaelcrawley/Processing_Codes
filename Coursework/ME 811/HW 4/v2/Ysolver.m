function [vh Resy] = Ysolver(alpha,AyO,AyE,AyW,AyN,AyS,Sy,v)
    Rtol = 1E-20;
    itrmax = 2;
    [M N] = size(Sy);
% %     Reshape link coefficients
%     AyO = reshape(AyO,N*M,1);
%     AyE = reshape(AyE,N*M,1);
%     AyW = reshape(AyW,N*M,1);
%     AyN = reshape(AyN,N*M,1);
%     AyS = reshape(AyS,N*M,1);
%     
%     AyN = circshift(AyN,-1);
%     AyS = circshift(AyS,1);%need to check to make sure AN,AS are being shifted properly
%     AyE = circshift(AyE,M);
%     AyW = circshift(AyW,-M);
%     
%     Ay = spdiags([AyW AyN AyO AyS AyE],[-M -1 0 1 M],N*M,N*M);
%     
% %     Modify source for correction equation
%     Ry = reshape(Sy,M*N,1)-Ay*reshape(v(1:end-1,:),M*N,1);
%     Sy = reshape(Ry,M,N);
%     Resy = norm(Ry);
%     Ay = Ay+diag(AyO*alpha);
%     
% %     Solve for v hat
%     [vh R] = ADI2d(Ay,Sy,'-TDMA',Rtol,itrmax);
    AyO = flipud([ones(1,N); AyO; ones(1,N)]);
    AyE = flipud([zeros(1,N); AyE; zeros(1,N)]);
    AyW = flipud([zeros(1,N); AyW; zeros(1,N)]);
    AyS = flipud([zeros(1,N); AyS; zeros(1,N)]);
    AyN = flipud([zeros(1,N); AyN; zeros(1,N)]);
    Sy = flipud([zeros(1,N); -Sy; zeros(1,N)]);
    v = flipud([zeros(1,N); v ]);
    
    [Resy,Ri] = calcRes(AyO,AyE,AyW,AyN,AyS,Sy,v);    
%     [vh,~] = ADIp((1+alpha)*AyO,AyE,AyW,AyN,AyS,Ri,Rtol,itrmax);
    [vh,~] = ADI(AyS,AyW,(1+alpha)*AyO,AyE,AyN,Ri,zeros(M+2,N),Rtol,itrmax);
    %append boundary condition
    vh = flipud(vh(2:end,:));
end