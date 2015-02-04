function out = Vj(M,varargin)
%Calculates the jet exit velocity (m/s)
%
%INPUTS
% M - mach number
% To - Stagnation total temperature (K)
% TTR - Total temperature ratio
% Ta - ambient temperature (K)

if nargin==2
    To = varargin{1};
    out = Tempj(M,To);
else
    TTR = varargin{1};
    Ta = varargin{2};
    out = Tempj(M,TTR,Ta);
end

out.Vj = M*sqrt(1.4*287.05*out.Tj);