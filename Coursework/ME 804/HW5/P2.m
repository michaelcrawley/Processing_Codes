%Problem 2, HW# 5, ME 804
%Michael Crawley
clear;clc;

T = 300;
Tv = 3000;
theta = 3395;
delta = 6e-3;
deltaVV = 6.5*T^(-1/2);
deltaVT = 2.1*T^(-1/3);
vend = 54;

%calc Boltzmann distribution
v = 0:vend;
fB = exp(-v*theta/Tv);

%calc Treaner distribution
vs = floor(0.5+T/2/Tv/delta);
v1 = 0:vs;
fT = exp(-theta*(1-2*delta)/Tv*v1+theta*delta/T*v1.*(v1-1));

%calc plateau-falloff
v2 = vs+1:vend;
k10 = 1.6e-20;
k0110 = 1e-14;
gamma = (vs+1)*(exp(-theta*(1-2*delta)/Tv*vs+theta*delta/T*vs*(vs-1))+(k10*T*(deltaVV^3)*exp(deltaVT*vs))/(k0110*12*theta*deltaVT*(vs+1)));
fP = gamma./(v2+1)-(k10*T*(deltaVV^3)*exp(deltaVT*v2))./(k0110*12*theta*deltaVT*(v2+1));

%combine
vt = [v1 v2];
fTP = [fT fP];

%calc local vibrational temperature
Tvl = theta*(1-2*delta*vt(2:end))./log(fTP(1:end-1)./fTP(2:end));

%plot
%code for dual-y-axis graph courtesy of Nathan Webb
h = figure; h1 = axes;
semilogy(v,fB,vt,fTP);
set(h1,'box','off');

h2 = axes;
plot(vt(2:end),Tvl);
set(h2,'box','off'); set(h2,'xaxislocation','top'); set(h2,'yaxislocation','right');
set(h2,'color','none'); set(h2,'xticklabel',{});




