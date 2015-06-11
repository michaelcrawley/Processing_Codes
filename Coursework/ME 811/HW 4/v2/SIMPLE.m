function [u v p Resx Resy] = SIMPLE(rho,gamma,alpha,omega,Ulid,M,N)
    Rtol = 1E-11;
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
    Resx = zeros(1,itrmax+1);Resx(1) = 1;
    counter = 1;
    Resy = Resx; Resp = Resx;
    
    while (Resx(counter) >= Rtol || Resy(counter) >= Rtol) && counter <= itrmax
        fprintf('Counter = %d\n',counter);
        fprintf('Calculating link coefs\n');
        [AxO AxE AxW AxN AxS Sx] = calcXcoefs(u,v,p,Ulid,rho,gamma,dx,dy);
        [AyO AyE AyW AyN AyS Sy] = calcYcoefs(u,v,p,rho,gamma,dx,dy);
        
        fprintf('Solving X momentum\n');
        [uh Resx(counter+1)] = Xsolver(alpha,AxO,AxE,AxW,AxN,AxS,Sx,u);
        uh = u+uh;
        fprintf('Solving Y momentum\n');
        [vh Resy(counter+1)] = Ysolver(alpha,AyO,AyE,AyW,AyN,AyS,Sy,v);
        vh = v+vh;
        
        fprintf('Calculating P link coefs\n');
        [ApO ApE ApW ApN ApS Sp] = calcPcoefs(rho,dx,dy,uh,vh,AxO,AyO);
        fprintf('Solving p correction\n');
        [pp Resp(counter+1)] = Psolver(ApO,ApE,ApW,ApN,ApS,Sp);
        
        fprintf('Updating values\n');
        u(:,2:end) = uh(:,2:end)+omega.u*dy*(pp(:,1:end-1)-pp(:,2:end))./AxO;
        v(1:end-1,:) = vh(1:end-1,:)+omega.v*dx*(pp(2:end,:)-pp(1:end-1,:))./AyO;
        p = p+omega.p*pp;
        
        counter = counter+1;
        clc;
    end
    u = [u zeros(M,1)];
    u = flipud(u);
    v = [zeros(1,N); v];
    v = flipud(v);
    p = flipud(p);
    Resx = Resx(2:counter);
    Resy = Resy(2:counter);
    h(1) = figure;contourf(x,yc,u);title('u velocity');colormap gray;colorbar;
    h(2) = figure;contourf(xc,y,v);title('v velocity');colormap gray;colorbar;
    h(3) = figure;contourf(xc,yc,p);title('pressure');colormap gray;colorbar;
    
    foldername = [num2str(M),'x',num2str(N),' Re',num2str(round(rho*Ulid*L/gamma))];
    if exist(foldername,'file') ~= 7
        mkdir(pwd,foldername);
    end
    cd(foldername);
    save data.mat;
    saveas(h(1),'u velocity','fig'); saveas(h(1),'u velocity','png');
    saveas(h(2),'v velocity','fig'); saveas(h(2),'v velocity','png');
    saveas(h(3),'p velocity','fig'); saveas(h(3),'p velocity','png');
    
    close(h);
    cd ..   
end