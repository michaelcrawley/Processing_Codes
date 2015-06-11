function [AxO AxE AxW AxN AxS Sx] = calcXcoefs(u,v,p,Ulid,rho,gamma,dx,dy)
%Calculates link coefficients for x-momentum FVM equation
%Note: for an NxMy mesh, u(:,1) is zero while v(:,1) and p(:,1) are not;
%the u and v velocity fields are modified to reflect this.
%Inputs:
%   u:      x-velocity component (NxM matrix)
%   v:      y-velocity component (NxM matrix)
%   p:      pressure (NxM matrix)
%   Ulid:   lid velocity
%   rho:    density (constant)
%   gamma:  viscosity (constant)
%   dx:     step size in x-direction
%   dy:     step size in y-direction


    [M N] = size(u);
    ut = [u zeros(M,1)]; %add known boundary values to simplify computations
    vt = [zeros(1,N); v];  %add known boundary values to simplify computations
    
    ce = rho*0.5*(ut(:,2:end-1)+ut(:,3:end));
    cw = rho*0.5*(ut(:,2:end-1)+ut(:,1:end-2));
    cn = rho*0.5*(vt(1:end-1,2:end)+vt(1:end-1,1:end-1));
    cs = rho*0.5*(vt(2:end,2:end)+vt(2:end,1:end-1));
    
    dep = 0.5*(abs(ce)+ce);
    dem = 0.5*(abs(ce)-ce);
    dwp = 0.5*(abs(cw)+cw);
    dwm = 0.5*(abs(cw)-cw);
    dsp = 0.5*(abs(cs)+cs);
    dsm = 0.5*(abs(cs)-cs);
    dnp = 0.5*(abs(cn)+cn);
    dnm = 0.5*(abs(cn)-cn);
    
    %Calculate link coefficients for interior nodes
    AxO = (dep+dwm+2*gamma/dx)*dy+(dnp+dsm+2*gamma/dy)*dx;
    AxE = -(dem+gamma/dx)*dy;
    AxW = -(dwp+gamma/dx)*dy;
    AxN = -(dnm+gamma/dy)*dx;
    AxS = -(dsp+gamma/dy)*dx;
    Sx = - (p(:,2:end)-p(:,1:end-1))*dy;
    
    %Modify link coefficients for North Boundary
    AxO(1,:) = (dep(1,:)+dwm(1,:)+2*gamma/dx)*dy+(dnp(1,:)+4*gamma/dy)*dx;
    AxN(1,:) = 0;
    AxS(1,:) = -(dsp(1,:)+4*gamma/dy/3)*dx;
    Sx(1,:) = Sx(1,:)+8*gamma*dx*Ulid/3/dy;
    
    %Modify link coefficients for South Boundary
    AxO(M,:) = (dep(M,:)+dwm(M,:)+2*gamma/dx)*dy+(dnp(M,:)+4*gamma/dy)*dx;
    AxS(M,:) = 0;
    AxN(M,:) = -(dnm(M,:)+4*gamma/dy/3)*dx;    
    
end