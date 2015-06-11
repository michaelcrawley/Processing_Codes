%ME 726 HW#4 Problem 12
%Written by Michael Crawley

clear;clc;

%define physical constants
R = 8.314; %in J/kmolK
Rb = 1.987; %in cal/molK
P = 101325; %in Pa
T1 = 500;
T2 = 3000;

%define equilibrium values
z1 = 5.5369e-9;
z2 = 0.3679;
n = 4+3*79/21;
yNOeq1 = z1/n;
yNOeq2 = z2/n;
yO2 = 1/n;
yN2 = 3*79/21/n;
C.NOeq1 = yNOeq1*P/R/T1;
C.NOeq2 = yNOeq2*P/R/T2;
C.Oeq1 = 7.158e-23;
C.Oeq2 = 0.1169;
C.O2eq1 = yO2*P/R/T1;
C.O2eq2 = yO2*P/R/T2;
C.N2eq1 = 9.6618e-046;
C.N2eq2 = 4.8407e-005;
C.Neq1 = 9.6618e-46;
C.Neq2 = 4.8407e-5;
C.Oeq1 = 7.158e-23;
C.Oeq2 = 0.1169;

%define reaction rate coefficients
k1_1 = (1.36e14)*exp(-75400/Rb/T1);
k1_2 = (1.36e14)*exp(-75400/Rb/T2);
k2_1 = (3.10e13)*exp(-334/Rb/T1);
k2_2 = (3.10e13)*exp(-334/Rb/T2);
k3_1 = (6.43e9)*T1*exp(-6250/Rb/T1);
k3_2 = (6.43e9)*T2*exp(-6250/Rb/T2);
k4_1 = k3_1*(6.5e-30);
k4_2 = k3_2*(1.04e-6);

%define equation coefs
a1 = k1_1*C.Oeq1*C.N2eq1+k3_1*C.Neq1*C.O2eq1;
a2 = k1_2*C.Oeq2*C.N2eq2+k3_2*C.Neq2*C.O2eq2;
b1 = -(k2_1*C.Neq1+k4_1*C.Oeq1);
b2 = -(k2_2*C.Neq2+k4_2*C.Oeq2);

t = 0:1e-6:1e-3;

NO1 = -(a2/b2)+(1/b2)*(a2+b2*C.NOeq1)*exp(b2*t);
NO2 = -(a1/b1)+(1/b1)*(a1+b1*C.NOeq2)*exp(b1*t);