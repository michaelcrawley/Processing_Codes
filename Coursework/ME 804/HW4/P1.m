%ME 804, Homework 4,Problem 1
%Written by Michael Crawley
clear;clc;

load thermo.mat

R = 8.314;
% k = 1.380e-23; %boltzman constant in J/K
k = 8.614e-5; %boltzman constant in eV/K
T = O.T;

%Calculate reaction coeff
dG = N.G+O.G-NOp.G-e.G;
krec = (4e-7)*(300^1.5)*T.^-1.5;
kion = krec.*exp(dG/R./T);

%Estimate Arrhenius rate

Arrhenius = @(x,E) ((4e-7)*300^1.5)*(x.^-1.5).*exp(-E./x/8.614e-5);
% E = GeneralFit(T,kion,Arrhenius);
E = 3.35;

plot(1./T,log(kion),'o',1./T,log(Arrhenius(T,E)));
title('Natural Logarithm of k_i_o_n vs 1/T');
xlabel('1/T (1/K)');
ylabel('ln k_i_o_n');
legend('k_i_o_n','Arrhenius rate expression','Location','NorthEast');