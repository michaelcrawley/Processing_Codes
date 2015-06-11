function [AyO AyE AyW AyN AyS Sy] = calcYcoefs(u,v,p,rho,gamma,dx,dy)
%Calculates link coefficients for y-momentum FVM equation
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
    
    ce = rho*0.5*(ut(1:end-1,2:end)+ut(2:end,2:end));
    cw = rho*0.5*(ut(1:end-1,1:end-1)+ut(2:end,1:end-1));
    cn = rho*0.5*(vt(2:end-1,:)+vt(1:end-2,:));
    cs = rho*0.5*( vt(2:end-1,:)+vt(3:end,:));
    
    dep = 0.5*(abs(ce)+ce);
    dem = 0.5*(abs(ce)-ce);
    dwp = 0.5*(abs(cw)+cw);
    dwm = 0.5*(abs(cw)-cw);
    dsp = 0.5*(abs(cs)+cs);
    dsm = 0.5*(abs(cs)-cs);
    dnp = 0.5*(abs(cn)+cn);
    dnm = 0.5*(abs(cn)-cn);
    
    %Calculate link coefficients for interior nodes
    AyO = (dep+dwm+2*gamma/dx)*dy+(dnp+dsm+2*gamma/dy)*dx;
    AyE = -(dem+gamma/dx)*dy;
    AyW = -(dwp+gamma/dx)*dy;
    AyN = -(dnm+gamma/dy)*dx;
    AyS = -(dsp+gamma/dy)*dx;
    Sy = -(p(1:end-1,:)-p(2:end,:))*dx;
    
    %Modify link coefficients for North Boundary
%     AyN(1,:) = 0;
    
    %Modify link coefficients for South Boundary - unnecessary
%     AyS(M-1,:) = 0;   
    
    %Modify link coefficients for East Boundary
    AyO(:,N) = (dwm(:,N)+4*gamma/dx)*dy+(dnp(:,N)+dsm(:,N)+2*gamma/dy)*dx;
    AyE(:,N) = 0; 
    AyW(:,N) = -(dwp(:,N)+4*gamma/dx/3)*dy;
    
    %Modify link coefficients for West Boundary
    AyO(:,1) = (dep(:,1)+4*gamma/dx)*dy+(dnp(:,1)+dsm(:,1)+2*gamma/dy)*dx;
    AyE(:,1) = -(dem(:,1)+4*gamma/dx/3)*dy;
    AyW(:,1) = 0;   
    
end