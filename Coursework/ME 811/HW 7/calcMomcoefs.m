function [Ax Ay] = calcMomcoefs(uf,vf,rho,P,gamma,Uin,dx,dy,xstepi,ystepi)
    [M N] = size(P);

    %Compute Face pressures
    p.e = [0.5*(P(:,1:end-1)+P(:,2:end)) zeros(M,1)]; %set outlet boundary condition
    p.w = [P(:,1) 0.5*(P(:,1:end-1)+P(:,2:end))];
    p.w(1:ystepi,xstepi+1) = P(1:ystepi,xstepi+1); %fix for step
    p.n = [0.5*(P(1:end-1,:)+P(2:end,:)); P(end,:)];
    p.s = [P(1,:); 0.5*(P(1:end-1,:)+P(2:end,:))];
    p.s(ystepi+1,1:xstepi) = P(ystepi+1,1:xstepi); %fix for step
    
    %%Calculate X-momentum link coefficients    
    ce = rho.e.*uf.e; %need to modify rho for weakly compressible flows
    cw = rho.w.*uf.w;
    cn = rho.n.*vf.n;
    cs = rho.s.*vf.s;
    
    dep = 0.5*(abs(ce)+ce);     
    dem = 0.5*(abs(ce)-ce);
    dwp = 0.5*(abs(cw)+cw);
    dwm = 0.5*(abs(cw)-cw);
    dnp = 0.5*(abs(cn)+cn);
    dnm = 0.5*(abs(cn)-cn);
    dsp = 0.5*(abs(cs)+cs);
    dsm = 0.5*(abs(cs)-cs);
    
    %Interior nodes
    Ax.O = (dep+dwm+2*gamma/dx)*dy+(dnp+dsm+2*gamma/dy)*dx;
    Ax.E = -(dem+gamma/dx)*dy;
    Ax.W = -(dwp+gamma/dx)*dy;
    Ax.N = -(dnm+gamma/dy)*dx;
    Ax.S = -(dsp+gamma/dy)*dx;
    Ax.P = -(p.e-p.w)*dy;
        
    %Inlet boundary
    Ax.O(ystepi+1:end,1) = Ax.O(ystepi+1:end,1)+2*gamma*dy/dx;
    Ax.E(ystepi+1:end,1) = Ax.E(ystepi+1:end,1)-gamma*dy/3/dx;
    Ax.W(ystepi+1:end,1) = 0;    
    Ax.P(ystepi+1:end,1) = Ax.P(ystepi+1:end,1)+rho.w(ystepi+1:end,1)*(Uin^2)*dy+8*gamma*Uin/3/dx*dy;
    
    %Ceiling boundary
    Ax.O(end,:) = Ax.O(end,:) - dnp(end,:)*dx+2*gamma*dx/dy;
    Ax.N(end,:) = 0;
    Ax.S(end,:) = Ax.S(end,:) - gamma*dx/3/dy;
    
    %Floor boundary
    Ax.O(1,:) = Ax.O(1,:)-dsm(1,:)*dx+2*gamma*dx/dy;
    Ax.N(1,:) = Ax.N(1,:) -gamma*dx/3/dy;
    Ax.S(1,:) = 0;
    
    %Outlet boundary
    Ax.O(:,end) = Ax.O(:,end) -dep(:,end)*dy-gamma*dy/dx+rho.e(:,end).*uf.e(:,end)*dy;
    Ax.E(:,end) = 0;
    
    %Step floor boundary
    Ax.O(ystepi+1,1:xstepi) = Ax.O(ystepi+1,1:xstepi)-dsm(ystepi+1,1:xstepi)*dx+2*gamma*dx/dy;
    Ax.N(ystepi+1,1:xstepi) = Ax.N(ystepi+1,1:xstepi) -gamma*dx/3/dy;
    Ax.S(ystepi+1,1:xstepi) = 0;
    
    %Step wall boundary
    Ax.O(1:ystepi,xstepi+1) = Ax.O(1:ystepi,xstepi+1)-dwm(1:ystepi,xstepi+1)*dy+2*gamma*dy/dx;
    Ax.E(1:ystepi,xstepi+1) = Ax.E(1:ystepi,xstepi+1)-gamma*dy/3/dx;
    Ax.W(1:ystepi,xstepi+1) = 0;
    
    %In step nodes
    Ax.O(1:ystepi,1:xstepi) = 1;
    Ax.E(1:ystepi,1:xstepi) = 0;
    Ax.W(1:ystepi,1:xstepi) = 0;
    Ax.N(1:ystepi,1:xstepi) = 0;
    Ax.S(1:ystepi,1:xstepi) = 0;
    Ax.P(1:ystepi,1:xstepi) = 0;
    
    %%Calculate Y-momentum link coefficients 
    Ay.O = Ax.O;
    Ay.E = Ax.E;
    Ay.W = Ax.W;
    Ay.S = Ax.S;
    Ay.N = Ax.N;
    Ay.P = -(p.n-p.s)*dx;
    %In step nodes
    Ay.P(1:ystepi,1:xstepi) = 0;
end