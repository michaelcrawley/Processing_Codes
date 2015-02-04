function alpha = atmAbsorption(f,T,P,RH)
% Function calculates sound attenuation coefficient by the atmosphere 
% according to ANSI S1.26(or ISO 9613-1) in dB/m. The ANSI and ISO
% standards are identical in their equations.

% List of input variables
%   f = frequency of sound wave (Hz)
%   T = measured temperature of ambient air (°C)
%   P = measured pressure of ambient air (kPa)
%   RH = measured relative humidity of ambient air (%)
% List of output variables
%   alpha = attenuation coefficient of air (dB/m)

T_ref = 293.15;       %Reference temperature (K)
P_ref = 101.325;      %Reference pressure (kPa)
T_01 = 273.16;        %Triple point isotherm temperature (K)

T_kel = T +273.15;     %Measured temperature conversion(K)
T_rel = T_kel/T_ref;  %Relative temperature
P_rel = P/P_ref;      %Relative pressure

%Saturation vapour pressure / reference pressure
V = 10.79586*(1-T_01/T_kel) -5.02808*log(T_kel/T_01) +1.50474e-4*(1-10^(-8.29692*(T_kel/T_01-1)))...
    -0.42873e-3*(1-10^(4.76955*(1-T_01/T_kel))) -2.2195983;
P_sat_over_P_ref = 10^V;

% % P_sat_over_P_ref = 10^(-6.8346*(T_01/T_kel)^1.261+4.6151);    %This is an approximate method which is accurate only near or below the triple point.

%h from ISO 9613-1, Annex B, eq.B1
H = RH*(P_sat_over_P_ref/P_rel);

%Oxygen relaxation frequency (Hz) from ISO 9613-1, 6.2, eq.3
fr_O = P_rel*(24+40400*H*(0.02+H)/(0.391+H));

%Nitrogen relaxation frequency (Hz) from ISO 9613-1, 6.2, eq.4
fr_N = P_rel/sqrt(T_rel)*(9+280*H*exp(-4.17*(T_rel^(-1/3)-1)));

%Intermediate calculations for classical, nitrogen and oxygen
%from ISO 9613-1, 6.2, part of eq.5
x_c = 1.84e-11/P_rel*sqrt(T_rel);
x_O = 0.01275*exp(-2239.1/T_kel)*fr_O./(fr_O^2+f.^2);
x_N = 0.1068*exp(-3352/T_kel)*fr_N./(fr_N^2+f.^2);

%Atmospheric attenuation coefficient, alpha (dB/m)
%from ISO 9613-1, 6.2, eq.5
alpha = 20*log10(exp(1))*f.^2.*(x_c +T_rel^(-5/2).*(x_O +x_N));
