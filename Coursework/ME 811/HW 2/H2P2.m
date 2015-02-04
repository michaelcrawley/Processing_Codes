function [phin Residual] = H2P2(N,M)
    x = 1/(2*N):1/N:1-1/(2*N);
    dx = mean(diff(x));
    y = 1/(2*M):1/M:1-1/(2*M);
    dy = mean(diff(y));
    n = 1:100; %number of terms for analytic solution
    [xx yy nn] = meshgrid(x,y,n);
    AR = dx/dy; iAR = dy/dx;
    
    %calculate analytic solution
    phia = (2/pi)*sum((((-1).^(nn+1)+1)./nn).*sin(nn*pi.*yy).*sinh(nn*pi.*xx)./sinh(nn*pi),3);
    
    gamma = ones(M,N,4); %last dimension is for the four faces; east, west, north, and south (in that order)
    S = zeros(M,N)*dx*dy;
    S(:,end) = S(:,end)-2*gamma(:,end,1)*iAR; %apply boundary conditions at x = 1
    gamma = reshape(gamma,M*N,4);
    
    X0 = -((gamma(:,1)+gamma(:,2))*iAR+(gamma(:,3)+gamma(:,4))*AR);
    X0(1:M) = X0(1:M)-gamma(1:M,2)*iAR; %apply boundary condition at x = 0
    X0(end-M+1:end) = X0(end-M+1:end)-gamma(end-M+1:end,1)*iAR; %apply boundary condition at x = 1
    X0(1:M:end) = X0(1:M:end)-gamma(1:M:end,4)*AR; %apply boundary condition at y = 0
    X0(M:M:end) = X0(M:M:end)-gamma(M:M:end,3)*AR; %apply boundary condition at y = 1
    X1 = gamma(:,3)*AR;
    X1(M:M:end) = 0; %apply boundary condition at y = 1
    Xn1 = gamma(:,4)*AR;
    Xn1(1:M:end) = 0; %apply boundary condition at y = 0
    XM = gamma(:,1)*iAR;
    XM(end-M+1:end) = 0; %apply boundary condition at x = 1
    XnM = gamma(:,2)*iAR;
    XnM(1:M) = 0; %apply boundary condition at x = 0;
    X1 = circshift(X1,1);
    XM = circshift(XM,M);
    Xn1 = circshift(Xn1,-1);
    XnM = circshift(XnM,-M);
    X = spdiags([XnM Xn1 X0 X1 XM],[-M -1 0 1 M],N*M,N*M);
    
    iResidual = norm(reshape(S,M*N,1));
    [phin Residual] = ADI2d(X,S,'-TDMA',iResidual/(10E6));
    figure; pcolor(x,y,phin); shading interp; colorbar; colormap(flipud(gray)); title(['Numeric Solution for ', num2str(N),'x',num2str(M),' nodes']);saveas(gcf,['Phin ', num2str(N),'x',num2str(M),' nodes'],'fig');
    figure; pcolor(x,y,phia-phin); shading interp; colorbar; colormap gray; title(['Numeric Error for ', num2str(N),'x',num2str(M),' nodes']);saveas(gcf,['Phie ', num2str(N),'x',num2str(M),' nodes'],'fig');
    figure; semilogy(Residual);title(['Residual for ', num2str(N),'x',num2str(M),' nodes']);xlabel('Iteration');ylabel('Residual');saveas(gcf,['Residual ', num2str(N),'x',num2str(M),' nodes'],'fig');

end