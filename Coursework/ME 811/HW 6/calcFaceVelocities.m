function [uf vf] = calcFaceVelocities(u,v,p,Uin,xstepi,ystepi,A,dx,dy,counter)
    [M N] = size(u);
    uf.e = zeros(M,N);
    uf.w = uf.e; vf.n = uf.e; vf.s = uf.e;

    if counter == 1 %DWIM
        uf.e = [0.5*(u(:,1:end-1)+u(:,2:end)) u(:,end)];
        uf.w = [Uin*ones(M,1) 0.5*(u(:,1:end-1)+u(:,2:end))];

        vf.n = [0.5*(v(1:end-1,:)+v(2:end,:)); zeros(1,N)];
        vf.s = [zeros(1,N); 0.5*(v(1:end-1,:)+v(2:end,:))];
    else %PWIM
        Ox = 3:N-2; Oy = 3:M-2;
        
        %interior nodes
        uf.e(:,Ox) = 0.5*(u(:,Ox)+u(:,Ox+1))+0.5*dx*dy*(1./A(:,Ox).*(p(:,Ox+1)-p(:,Ox-1))/2/dx+1./A(:,Ox+1).*(p(:,Ox+2)-p(:,Ox))/2/dx-(1./A(:,Ox)+1./A(:,Ox+1)).*(p(:,Ox+1)-p(:,Ox))/dx);
        uf.w(:,Ox) = 0.5*(u(:,Ox)+u(:,Ox-1))+0.5*dx*dy*(1./A(:,Ox).*(p(:,Ox+1)-p(:,Ox-1))/2/dx+1./A(:,Ox-1).*(p(:,Ox)-p(:,Ox-2))/2/dx-(1./A(:,Ox)+1./A(:,Ox-1)).*(p(:,Ox)-p(:,Ox-1))/dx);
        
        vf.n(Oy,:) = 0.5*(v(Oy,:)+v(Oy+1,:))+0.5*dx*dy*(1./A(Oy,:).*(p(Oy+1,:)-p(Oy-1,:))/2/dy+1./A(Oy+1,:).*(p(Oy+2,:)-p(Oy,:))/2/dy-(1./A(Oy,:)+1./A(Oy+1,:)).*(p(Oy+1,:)-p(Oy,:))/dy);        
        vf.s(Oy,:) = 0.5*(v(Oy,:)+v(Oy-1,:))+0.5*dx*dy*(1./A(Oy,:).*(p(Oy+1,:)-p(Oy-1,:))/2/dy+1./A(Oy-1,:).*(p(Oy,:)-p(Oy-2,:))/2/dy-(1./A(Oy,:)+1./A(Oy-1,:)).*(p(Oy,:)-p(Oy-1,:))/dy);
        
        %Near Outlet nodes
        uf.w(:,N-1) = uf.e(:,N-2);
        uf.e(:,N-1) = 0.5*(u(:,N-1)+u(:,N))+0.5*dx*dy*(1./A(:,N-1).*(p(:,N)-p(:,N-2))/2/dx+1./A(:,N).*(p(:,N)-p(:,N-1))/dx-(1./A(:,N-1)+1./A(:,N)).*(p(:,N)-p(:,N-1))/dx);
        
        %Outlet nodes
        uf.e(:,N) = u(:,N);
        uf.w(:,N) = uf.e(:,N-1);
        
        %Near inlet nodes
        uf.e(ystepi+1:end,2) = uf.w(ystepi+1:end,3);
        uf.w(ystepi+1:end,2) = 0.5*(u(ystepi+1:end,2)+u(ystepi+1:end,1))+0.5*dx*dy*(1./A(ystepi+1:end,2).*(p(ystepi+1:end,2+1)-p(ystepi+1:end,2-1))/2/dx+1./A(ystepi+1:end,2-1).*(p(ystepi+1:end,2)-p(ystepi+1:end,1))/dx-(1./A(ystepi+1:end,2)+1./A(ystepi+1:end,2-1)).*(p(ystepi+1:end,2)-p(ystepi+1:end,2-1))/dx);
        
        %Inlet nodes
        uf.w(ystepi+1:end,1) = Uin;
        uf.e(ystepi+1:end,1) = uf.w(ystepi+1:end,2);
        
        %Near step wall nodes
        uf.e(1:ystepi,xstepi+2) = uf.w(1:ystepi,xstepi+3);
        uf.w(1:ystepi,xstepi+2) = 0.5*(u(1:ystepi,xstepi+2)+u(1:ystepi,xstepi+2-1))+0.5*dx*dy*(1./A(1:ystepi,xstepi+2).*(p(1:ystepi,xstepi+2+1)-p(1:ystepi,xstepi+2-1))/2/dx+1./A(1:ystepi,xstepi+2-1).*(p(1:ystepi,xstepi+2)-p(1:ystepi,xstepi+2-1))/dx-(1./A(1:ystepi,xstepi+2)+1./A(1:ystepi,xstepi+2-1)).*(p(1:ystepi,xstepi+2)-p(1:ystepi,xstepi+2-1))/dx);
        
        %Step wall nodes
        uf.w(1:ystepi,xstepi+1) = 0;
        uf.e(1:ystepi,xstepi+1) = uf.w(1:ystepi,xstepi+2);
        
        %Near ceiling nodes
        vf.s(end-1,:) = vf.n(end-2,:);
        vf.n(end-1,:) = 0.5*(v(end-1,:)+v(end-1+1,:))+0.5*dx*dy*(1./A(end-1,:).*(p(end-1+1,:)-p(end-1-1,:))/2/dy+1./A(end-1+1,:).*(p(end-1+1,:)-p(end-1,:))/dy-(1./A(end-1,:)+1./A(end-1+1,:)).*(p(end-1+1,:)-p(end-1,:))/dy);
        
        %Ceiling nodes
        vf.n(end,:) = 0;
        vf.s(end,:) = vf.n(end-1,:);
        
        %Near step floor nodes
        vf.n(ystepi+2,1:xstepi) = vf.s(ystepi+3,1:xstepi);
        vf.s(ystepi+2,1:xstepi) = 0.5*(v(ystepi+2,1:xstepi)+v(ystepi+2-1,1:xstepi))+0.5*dx*dy*(1./A(ystepi+2,1:xstepi).*(p(ystepi+2+1,1:xstepi)-p(ystepi+2-1,1:xstepi))/2/dy+1./A(ystepi+2-1,1:xstepi).*(p(ystepi+2,1:xstepi)-p(ystepi+2-1,1:xstepi))/dy-(1./A(ystepi+2,1:xstepi)+1./A(ystepi+2-1,1:xstepi)).*(p(ystepi+2,1:xstepi)-p(ystepi+2-1,1:xstepi))/dy);
        
        %Step floor nodes
        vf.s(ystepi+1,1:xstepi) = 0;
        vf.n(ystepi+1,1:xstepi) = vf.s(ystepi+2,1:xstepi);
        
        %Near floor nodes
        vf.n(2,xstepi+1:end) = vf.s(3,xstepi+1:end);
        vf.s(2,xstepi+1:end) = 0.5*(v(2,xstepi+1:end)+v(2-1,xstepi+1:end))+0.5*dx*dy*(1./A(2,xstepi+1:end).*(p(2+1,xstepi+1:end)-p(2-1,xstepi+1:end))/2/dy+1./A(2-1,xstepi+1:end).*(p(2,xstepi+1:end)-p(2-1,xstepi+1:end))/dy-(1./A(2,xstepi+1:end)+1./A(2-1,xstepi+1:end)).*(p(2,xstepi+1:end)-p(2-1,xstepi+1:end))/dy);
        
        %Floor nodes
        vf.s(1,xstepi+1:end) = 0;
        vf.n(1,xstepi+1:end) = vf.s(2,xstepi+1:end);        
    end
    
    %modifications for step
    uf.w(1:ystepi,1:xstepi+1) = 0;
    uf.e(1:ystepi,1:xstepi) = 0;
    
    vf.s(1:ystepi+1,1:xstepi) = 0;
    vf.n(1:ystepi,1:xstepi) = 0;
end