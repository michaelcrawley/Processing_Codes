%Solution for ME 804 HW#2, Problem 2.  Written by Michael Crawley

clear;clc;

%Load thermo data from mat files
load 'Othermo.mat';
load 'O2thermo.mat';

%initialize constants
R = 8.314;
L = length(O.T);

%Part A
Kp = exp(-(2*O.G-O2.G)./O.T/R);
h(1) = figure;
plot(O.T,log10(Kp)); xlabel('T (K)'); ylabel('Log_1_0 K_p'); title('Equilibrium constant for O_2 = 2O versus Temperature');
saveas(h(1),'Kp','fig');

%Part B
%calculate mole fractions
xO = 0.5*(-Kp+sqrt(Kp.^2+4*Kp));
xO2 = 1-xO;
h(2) = figure;
plot(O.T,xO); xlabel('T (K)'); ylabel('x_O'); title('Mole fraction of O versus temperature');
saveas(h(2),'mole fraction','fig');
%calculate degree of dissociation
alpha = zeros(L,1);
for i = 1:L
    n = null([xO(i)-1 xO(i); xO2(i) xO2(i)-1]);
    alpha(i) = n(1)/(n(1)+2*n(2));
end
h(3) = figure;
plot(O.T,alpha); xlabel('T (K)'); ylabel('\alpha'); title('Degree of Dissociation of O versus temperature');
saveas(h(3),'alpha','fig');

%Part C
H = xO.*O.H+xO2.*O2.H;
h(4) = figure;
plot(O.T,H); xlabel('T (K)'); ylabel('H (kJ/kmol)'); title('Enthalpy of mixture versus temperature');
saveas(h(4),'enthalpy','fig');

%Part D
Cpf = xO.*O.Cp+xO2.*O2.Cp;
h(5) = figure;
plot(O.T,Cpf); xlabel('T (K)'); ylabel('C_p_,_f_r_o_z_e_n'); title('C_p_,_f_r_o_z_e_n versus temperature');
saveas(h(5),'Cpf','fig');

%Part E
ddT = finite_diff_Taylor_non_uniform(O.T,1,2);
dxO = ddT*xO;
dxO2 = ddT*xO2;
Cpe = Cpf+O.H.*dxO+O2.H.*dxO2;
h(6) = figure;
plot(O.T,Cpe); xlabel('T (K)'); ylabel('C_p_,_e_q_u_i_l'); title('C_p_,_e_q_u_i_l versus temperature');
saveas(h(6),'Cpe','fig');

close all;