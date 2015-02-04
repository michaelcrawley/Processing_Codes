function [u v p Resx Resy Resp] = SIMPLEcol(Uin,M)
    tic;
    Rtol = 10E-10;
    itrmax = 2E3;
    
    N = 5*M; %force uniform dx, dy
    rho = ones(M,N);
    gamma = 2E-5;
    H = 0.01; %m
    
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
    p = v; p(:,end) = pout;

    Resx = zeros(1,itrmax+1);Resx(1) = 1;
    counter = 1;
    Resy = Resx; Resp = Resx;
    Ax = []; 
    [uf vf] = calcFaceVelocities(u,v,p,Uin,xstepi,ystepi,Ax,dx,dy,counter);
    while (Resx(counter) >= Rtol || Resy(counter) >= Rtol || Resp(counter) >= Rtol) && counter <= itrmax
        counter = counter+1;
        [Ax Ay] = calcMomcoefs(uf,vf,rho,p,gamma,Uin,dx,dy,xstepi,ystepi);
        [u Resx(counter)] = solveMomEQ(u,Ax,alpha);
        [v Resy(counter)] = solveMomEQ(v,Ay,alpha);
        
        [uf vf] = calcFaceVelocities(u,v,p,Uin,xstepi,ystepi,Ax.O,dx,dy,counter);
        [Ap] = calcPresscoefs(uf,vf,rho,dx,dy,xstepi,ystepi,Ax.O);
        [pp Resp(counter)] = solvePressEQ(Ap);        
        
        [u,uf,v,vf,p] = updatevalues(u,uf,v,vf,p,pp,dx,dy,xstepi,ystepi,omega,Ax.O);        
        
        fprintf('counter: %i Resx: %1.2e Resy: %1.2e Resp: %1.2e \n',[counter-1 Resx(counter) Resy(counter) Resp(counter)]);                
    end
    
    Re=round(mean(mean(rho))*Uin*2*H/gamma);
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
    
    %Save figures
    saveas(h(1),[filename,' pressure'],'fig');saveas(h(1),[filename,' pressure'],'png');
    saveas(h(2),[filename,' Residuals'],'fig');saveas(h(2),[filename,' Residuals'],'png');
    saveas(h(3),[filename,' u velocity'],'fig');saveas(h(3),[filename,' u velocity'],'png');    
    saveas(h(4),[filename,' v velocity'],'fig');saveas(h(4),[filename,' v velocity'],'png');

    close(h);
    cd .. 
end