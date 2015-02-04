function [At] = calcTcoefs(uf,vf,rhof,p,Uin,dx,dy,xstepi,ystepi,stepT,cp,k)
    rhof.e = rhof.e*cp;
    rhof.w = rhof.w*cp;
    rhof.n = rhof.n*cp;
    rhof.s = rhof.s*cp;
    
    [At,~] = calcMomcoefs(uf,vf,rhof,p,k,Uin,dx,dy,xstepi,ystepi);
    
    At.P = zeros(size(At.P));
    At.P(ystepi+1,1:xstepi) = 8*k*stepT/3;
    At.P(1:ystepi,xstepi+1) = 8*k*stepT/3;
end