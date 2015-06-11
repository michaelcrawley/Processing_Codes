%ME 804, Homework 4,Problem 2
%Written by Michael Crawley
clear;clc

load O.mat
load O2.mat

R = 8.314;
Na = 6.022E23; %molecule/mole
k = 1.380e-23;

P = 30.75*101325; %pressure, in Pa
T = 300*8.3; %temperature, in K
M = 0.468;
U = M*sqrt(1.63*R*1000/39.2*T);
I = find(abs(O.T-T) == min(abs(O.T-T)));
C_O2 = 0.1*P/R/T;
C_Ar = 0.9*P/R/T;

%Arrhenius rate coeffs
A = (1e-8)*(100^-3)*Na; %m3/mole/s
Ea = 494000; %in J

kf = A*(298/T)*exp(-Ea/R/T);
keq = exp(-(2*O.G(I)-O2.G(I))/R/T)*(P/101325);
kr = kf/keq;


t = 0:1e-6:2e-2;
s = U*t;
xO = (R*T/P)*sqrt(keq*C_O2)*tanh(C_Ar*t*sqrt(kf*kr*C_O2));
xOeq = max(xO);
I = find(xO > 0.99*xOeq,1,'first');

plot(s,xO,s(I),xO(I),'o');
title('O mole fraction vs distance behind shock');
xlabel('m');
ylabel('y_o');
