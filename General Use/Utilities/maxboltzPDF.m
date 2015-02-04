function y = maxboltzPDF(x,S,varargin)
%Function creates the Maxwell-Boltzmann probability density function one of
%two ways. If the only input arguments are the x-axis (x) and the width
%parameter (S), it produces the traditional PDF which integrates to 1. If
%there are more arguments, it assumes the following inputs:
%
%  x: x-axis
%  S: width - often written as sigma in symbolic representations
% xo: horizontal offset
% yo: vertical offset
%  A: amplitude
%
%The program uses a peak value normalized form of the PDF so that an
%arbitray amplitude may be applied as well as allowing for an arbitrary
%shift of origin. In this case, a heaviside function is also used to
%eliminate y-axis symmetry.

if isempty(varargin)
    y = sqrt(2/pi)*x.^2.*exp(-x.^2/2/S^2)/S^3;
else
    xo = varargin{1};
    yo = varargin{2};
    A = varargin{3};
    
    y = A*heaviside(x-xo).*(x-xo).^2/2/S^2.*exp(1-(x-xo).^2/2/S^2) +yo;
end
