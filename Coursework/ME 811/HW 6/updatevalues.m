function [u,uf,v,vf,p] = updatevalues(u,uf,v,vf,p,pp,dx,dy,xstepi,ystepi,omega,Ao)
    [M N] = size(u);

    %update cell center pressure
    p = p+omega.p*pp;
    
    ppx = [pp(:,1) pp zeros(M,1)];
    ppy = [pp(1,:); pp; pp(end,:)];
    
    ppx(1:ystepi,xstepi+1) = ppx(1:ystepi,xstepi+2);
    ppy(ystepi+1,1:xstepi) = ppy(ystepi+2,1:xstepi);
        
    %update cell center velocities
    up = (ppx(:,1:end-2)-ppx(:,3:end))*dy/2./Ao; 
    vp = (ppy(1:end-2,:)-ppy(3:end,:))*dx/2./Ao;
    
    up(1:ystepi,1:xstepi) = 0;
    vp(1:ystepi,1:xstepi) = 0;
    
    u = u+omega.u*up;
    v = v+omega.u*vp;
    
    %update cell face velocities
    uf.ep = zeros(size(p));
    uf.wp = uf.ep;
    vf.np = uf.ep;
    vf.sp = uf.ep;
    
    uf.ep(:,1:end-1) = 0.5*dy*(1./Ao(:,1:end-1)+1./Ao(:,2:end)).*(pp(:,1:end-1)-pp(:,2:end));
    uf.wp(:,2:end) = 0.5*dy*(1./Ao(:,2:end)+1./Ao(:,1:end-1)).*(pp(:,1:end-1)-pp(:,2:end));
    
    uf.ep(1:ystepi,xstepi) = 0;
    uf.wp(1:ystepi,xstepi+1) = 0;
    
    vf.np(1:end-1,:) = 0.5*dx*(1./Ao(1:end-1,:)+1./Ao(2:end,:)).*(pp(1:end-1,:)-pp(2:end,:));
    vf.sp(2:end,:) = 0.5*dx*(1./Ao(2:end,:)+1./Ao(1:end-1,:)).*(pp(1:end-1,:)-pp(2:end,:));
    
    vf.np(ystepi,1:xstepi) = 0;
    vf.sp(ystepi+1,1:xstepi) = 0;
    
    uf.e = uf.e+omega.u*uf.ep;
    uf.w = uf.w+omega.u*uf.wp;
    vf.n = vf.n+omega.u*vf.np;
    vf.s = vf.s+omega.u*vf.sp;
end