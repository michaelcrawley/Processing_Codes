function [rhoc] = calcRho(pc,Tc,pref,Tref,R)

[M N] = size(pc);

i=1:M; j=1:N;
rhoc(i,j) = 1/R * (pref+pc(i,j))./(Tref+Tc(i,j));