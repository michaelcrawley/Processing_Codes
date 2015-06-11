function out = Tempj(M,varargin)
%Calculates the jet exit Temp (K)
%
%INPUTS
% M - mach number
% To - Stagnation total temperature (K)
% TTR - Total temperature ratio
% Ta - ambient temperature (K)

if nargin==2
    To = varargin{1};
    TTR = [];
    Ta = [];
else
    TTR = varargin{1};
    Ta = varargin{2};
    To = TTR*Ta;
end

Tj = To/(1+0.2*M^2);

out.M = M;
out.TTR = TTR;
out.Ta = Ta;
out.To = To;
out.Tj = Tj;