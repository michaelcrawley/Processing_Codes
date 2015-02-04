function [rho rhof] = calcrho(T,p,Tref,pref,R)

    rho = (1/R)*(p+pref)./(T+Tref);
    rhof.e = 0.5*(rho+[rho(:,2:end) rho(:,end)]);
    rhof.w = 0.5*(rho+[rho(:,1) rho(:,1:end-1)]);
    rhof.n = 0.5*(rho+[rho(2:end,:); rho(end,:)]);
    rhof.s = 0.5*(rho+[rho(1,:); rho(1:end-1,:)]);
    
end