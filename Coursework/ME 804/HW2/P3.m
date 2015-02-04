%Solution for ME 804 HW#2, Problem 3.  Written by Michael Crawley

clear;clc;

I = 3.89*1.602e-19; %ionization energy of cesium, in J
k = 1.380e-23; %boltman constant, in J/K
h = 6.626e-34; %planck constant, in J*s
R = 8.314; %gas constant, in J/mol/K
P = 0.01*101325; %pressure, in Pa
me = 9.109e-31; %electron mass, in kg
mcs = 132.9*1.661e-27; %cesium atom mass, in kg
M = 132.9; %cesium molar mass

%degeneracies
gcs = 2;
gcsp = 1;
ge = 2;

%define temperature range
T = 300:1:6000;

rho = P*M/R./T;
RHS = (mcs*gcsp*ge/gcs)*(1./rho).*((2*pi*me*k*T/h/h).^(3/2)).*exp(-I/k./T);
phi = 0.5*(-RHS+sqrt(RHS.^2+4*RHS)); %calculate ionization fraction
xe = phi./(1+phi);

plot(T,phi,'k',T,xe,'--k');legend('\phi','x_e','Location','NorthWest'); xlabel('Temperature (K)'); ylabel('Level');title({'Degree of Ionization and mole fraction of electrons','versus temperature'});
saveas(gcf,'P2','fig');close all;