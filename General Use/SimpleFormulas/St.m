function out = St(f,L,U)
%Calculates the Strouhal number.
%
%INPUTS
% f - frequency (Hz)
% L - length scale (m)
% U - velocity scale (m/s)

out = f*L/U;