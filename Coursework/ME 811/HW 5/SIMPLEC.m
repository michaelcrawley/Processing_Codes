function [u v p Resx Resy] = SIMPLEC(rho,gamma,alpha,omega,Ulid,M,N)
    tic
    Rtol = 10E-10;
    itrmax = 1E4;
    L = 0.01;
    dx = L/N;
    dy = L/M;
    x = 0:dx:L;
    y = 0:dy:L;
    xc = dx/2:dx:L-dx/2;
    yc = dy/2:dy:L-dy/2;

    %initialize variables
    u = zeros(M,N);
    v = u; p = u;
%     v = ones(M,N); p = v; v(end,:) = 0; u = ones(M,N); u(:,1) = 0;
    Resx = zeros(1,itrmax+1);Resx(1) = 1;
    counter = 1;
    Resy = Resx; Resp = Resx;
    
    while (Resx(counter) >= Rtol || Resy(counter) >= Rtol || Resp(counter) >= Rtol) && counter <= itrmax
        [AxO AxE AxW AxN AxS Sx] = calcXcoefs(u,v,p,Ulid,rho,gamma,dx,dy);
        [AyO AyE AyW AyN AyS Sy] = calcYcoefs(u,v,p,rho,gamma,dx,dy);
        
        [uh Resx(counter+1)] = Xsolver(alpha,AxO,AxE,AxW,AxN,AxS,Sx,u);
        uh = u+uh;

        [vh Resy(counter+1)] = Ysolver(alpha,AyO,AyE,AyW,AyN,AyS,Sy,v);
        vh = v+vh;
        
        [ApO ApE ApW ApN ApS Sp] = calcPcoefs(alpha,rho,dx,dy,uh,vh,AxO,AxE,AxW,AxN,AxS,AyO,AyE,AyW,AyN,AyS);
        [pp Resp(counter+1)] = Psolver(ApO,ApE,ApW,ApN,ApS,Sp);
        
        u(:,2:end) = uh(:,2:end)+omega*dy*(pp(:,1:end-1)-pp(:,2:end))./((1+alpha)*AxO+AxE+AxW+AxN+AxS);
        v(1:end-1,:) = vh(1:end-1,:)+omega*dx*(pp(2:end,:)-pp(1:end-1,:))./((1+alpha)*AyO+AyE+AyW+AyN+AyS);
        p = p+pp;
        
        fprintf('counter: %i Resx: %1.2e Resy: %1.2e Resp: %1.2e \n',[counter Resx(counter+1) Resy(counter+1) Resp(counter+1)]);
        counter = counter+1;
    end
    u = [u zeros(M,1)];
    u = flipud(u);
    v = [zeros(1,N); v];
    v = flipud(v);
    p = flipud(p);
    Resx = Resx(2:counter);
    Resy = Resy(2:counter);
    Resp = Resp(2:counter);
    
    foldername = [num2str(M),'x',num2str(N),' Re',num2str(round(rho*Ulid*L/gamma))];
    if exist(foldername,'file') ~= 7
        mkdir(pwd,foldername);
    end
    cd(foldername);
    compute_time = toc;
    save data.mat;

    filename = [num2str(M),'x',num2str(N),' mesh, Re = ',num2str(round(rho*Ulid*L/gamma))];

    %Plot centerline u,v velocity along centerlines
    h(1) = figure; %y = const figure
    plot(x,(u(M/2,:)+u(M/2+1,:))/2,'k',xc,v(M/2+1,:),'k--');legend('u velocity','v velocity','Location','Best');xlabel('x');ylabel('(m/s)');title(['Horizontal Centerline Velocities for ',filename]);
    h(2) = figure; %x = const figure
    plot(yc,u(:,N/2+1),'k',y,(v(:,N/2)+v(:,N/2+1))/2,'k--');legend('u velocity','v velocity','Location','Best');xlabel('y');ylabel('(m/s)');title(['Vertical Centerline Velocities for ',filename]);

    %Plot residuals
    h(3) = figure;
    semilogy(1:counter-1,Resx,'k',1:counter-1,Resy,'k--',1:counter-1,Resp,'k-.');legend('u Residual','v Residual','p Residual');xlabel('iteration');ylabel('Residual');title(['Outer Iteration Residuals for ',filename]);

    %Plot quiver and pressure
    uc = (u(:,1:end-1)+u(:,2:end))/2;
    vc = (v(1:end-1,:)+v(2:end,:))/2;
    h(4) = figure;
    contourf(xc,yc,p,50);colormap(flipud(gray));colorbar;hold on;
    quiver(xc,yc,uc,vc,3);xlabel('x');ylabel('y');hold off;
    title(['Velocity Vectors and Pressure Contours for ',filename]);

    %Save figures
    saveas(h(1),[filename,' y centerline'],'fig');saveas(h(1),[filename,' y centerline'],'png');
    saveas(h(2),[filename,' x centerline'],'fig');saveas(h(2),[filename,' x centerline'],'png');
    saveas(h(3),[filename,' Residuals'],'fig');saveas(h(3),[filename,' Residuals'],'png');
    saveas(h(4),[filename,' pressure'],'fig');saveas(h(4),[filename,' pressure'],'png');

    close(h);
    cd ..  
end