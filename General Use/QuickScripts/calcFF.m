function [ForcingFrequency]=calcFF(t0, M, StD,X)
%Inputs: To, M, St, X
T0=t0+273;%Convert to degrees kelvin
Te=T0/(1+0.2*M^2);
c=(1.4*287*Te)^0.5;
U=M*c;
ForcingFrequency=(StD*U/X);
end
