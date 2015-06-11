function out = StInv(St,L,U)
%Calculates the frequency from the Strouhal number.
%
%INPUTS
% f - frequency (Hz)
% L - length scale (m)
% U - velocity scale (m/s)

out = St*U/L;