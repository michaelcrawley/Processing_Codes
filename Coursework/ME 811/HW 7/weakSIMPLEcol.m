function [u v p Resx Resy Resp] = weakSIMPLEcol(M)
    tic;
    Rtol = 10E-10;
    itrmax = 2E3;
    
    N = 5*M; %force uniform dx, dy
    mu = 2E-5;
    k = 0.026;
    cp = 1012;
    R = 287;
    Tref = 300;
    pref = 101325;
    Uin = 0.1;
    H = 0.01; %m
    stepT = 100;
    
    omega.u = 0.5;						
    omega.p = 0.2;						
    alpha = 0.2;	
    
    dx = 10*H/N;
    dy = 2*H/M;
    
    x = dx/2:dx:10*H-dx/2;
    y = dy/2:dy:2*H-dy/2;
    
    xstepi = find(x <= 2*H,1,'last'); %locate x-indices inside step
    ystepi = find(y <= H,1,'last'); %lcoate y-indices inside step
    
    pout = 0;
    u = Uin*ones(M,N); u(1:ystepi,1:xstepi) = 0;
    v = zeros(M,N); 
    T = v;
    p = v; p(:,end) = pout;

    Resx = zeros(1,itrmax+1);Resx(1) = 1;
    counter = 1;
    Resy = Resx; Resp = Resx; ResT = Resx;
    Ax = []; 
    [uf vf] = calcFaceVelocities(u,v,p,Uin,xstepi,ystepi,Ax,dx,dy,counter);
    [~, rhof] = calcrho(T,p,Tref,pref,R);
    while (Resx(counter) >= Rtol || Resy(counter) >= Rtol || Resp(counter) >= Rtol || ResT(counter) >= Rtol) && counter <= itrmax
        counter = counter+1;
        [Ax Ay] = calcMomcoefs(uf,vf,rhof,p,mu,Uin,dx,dy,xstepi,ystepi);
        [u Resx(counter)] = solveMomEQ(u,Ax,alpha);
        [v Resy(counter)] = solveMomEQ(v,Ay,alpha);
        
        [uf vf] = calcFaceVelocities(u,v,p,Uin,xstepi,ystepi,Ax.O,dx,dy,counter);
        [Ap] = calcPresscoefs(uf,vf,rhof,dx,dy,xstepi,ystepi,Ax.O);
        [pp Resp(counter)] = solvePressEQ(Ap);        
        
        [u,uf,v,vf,p] = updatevalues(u,uf,v,vf,p,pp,dx,dy,xstepi,ystepi,omega,Ax.O);
        
        [At] = calcTcoefs(uf,vf,rhof,p,Uin,dx,dy,xstepi,ystepi,stepT,cp,k);
        
        [T ResT(counter)]  = solveTEQ(T,At,alpha);
        
        [~, rhof] = calcrho(T,p,Tref,pref,R);
        
        fprintf('counter: %i Resx: %1.2e Resy: %1.2e Resp: %1.2e RespT: %1.2e\n',[counter-1 Resx(counter) Resy(counter) Resp(counter) ResT(counter)]);                
    end
    
    T = T+Tref; T(1:ystepi,1:xstepi) = T(1:ystepi,1:xstepi)+100;
    
    Re=round(1*Uin*2*H/mu);
    foldername = [num2str(M),'x',num2str(N),' Re',num2str(Re)];
    filename = [num2str(M),'x',num2str(N),' mesh, Re = ',num2str(Re)];
    if exist(foldername,'file') ~= 7
        mkdir(pwd,foldername);
    end
    cd(foldername);
    compute_time = toc;
    save data.mat;
    
    h(1) = figure;
    contourf(x,y,p,50);colormap(flipud(gray));colorbar;hold on;
    quiver(x,y,u,v,3);xlabel('x');ylabel('y');xlim([x(1) x(end)]); ylim([y(1) y(end)]);plot([0 0.02 0.02],[0.01 0.01 0],'k');hold off;
    xlabel('x (m)'); ylabel('y (m)'); title(['Velocity Vectors and Pressure Contours for ',filename]);
    
    h(2) = figure;
    semilogy(1:counter,Resx(1:counter),'k',1:counter,Resy(1:counter),'k--',1:counter,Resp(1:counter),'k-.');
    legend('u Residual','v Residual','p Residual');xlabel('iteration');ylabel('Residual');title(['Outer Iteration Residuals for ',filename]);
    
    h(3) = figure;
    contourf(x,y,u);hold on;plot([0 0.02 0.02],[0.01 0.01 0],'k');hold off;colormap(flipud(gray));colorbar;xlabel('x');ylabel('y');title(['U velocity Contours for ',filename]);
    
    h(4) = figure;
    contourf(x,y,v);hold on;plot([0 0.02 0.02],[0.01 0.01 0],'k');hold off;colormap(flipud(gray));colorbar;xlabel('x');ylabel('y');title(['V velocity Contours for ',filename]); 
    
    h(5) = figure;
    contourf(x,y,T);hold on;plot([0 0.02 0.02],[0.01 0.01 0],'k');hold off;colormap(flipud(gray));colorbar;xlabel('x');ylabel('y');title(['Temperature Contours for ',filename]); 
    
    %Save figures
    saveas(h(1),[filename,' pressure'],'fig');saveas(h(1),[filename,' pressure'],'png');
    saveas(h(2),[filename,' Residuals'],'fig');saveas(h(2),[filename,' Residuals'],'png');
    saveas(h(3),[filename,' u velocity'],'fig');saveas(h(3),[filename,' u velocity'],'png');    
    saveas(h(4),[filename,' v velocity'],'fig');saveas(h(4),[filename,' v velocity'],'png');
    saveas(h(5),[filename,' temperature'],'fig');saveas(h(5),[filename,' temperature'],'png');

    close(h);
    cd .. 
end