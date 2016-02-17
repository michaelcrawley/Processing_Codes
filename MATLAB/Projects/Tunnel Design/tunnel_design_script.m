
clear;

%Variables
M = 2.5;
A_test = (6*0.0254)^2;
V = 8000 * 0.00378541; % convert from gallons to cubic meters
p_i = 190*6894.76;

%Constants
R = 287;
n = 1.2;
gamma = 1.4;
p_r = 1.25; %pressure ratio for pressure drop
T = 300;
p_amb = 101325;

%Equations
gamma_mod = 0.5*(gamma+1)/(gamma-1);
A_r = ((gamma+1)/2)^(-gamma_mod) * ((1+ 0.5*(gamma-1)*M^2)^(gamma_mod))/M;
A_s = A_test/A_r;
p_1 = p_amb*(gamma+1)/(2*gamma*M^2 - (gamma-1));
p_t = p_1 * (1 + 0.5*(gamma-1)*M^2)^(gamma/(gamma-1));
p_f = p_r*p_t;
t = ((1 + (gamma-1)/2)^gamma_mod) * sqrt(1/gamma/R) * (1/sqrt(T)) * (p_i/p_t) * (V/A_s) * (1 - (p_f/p_i)^(1/n));