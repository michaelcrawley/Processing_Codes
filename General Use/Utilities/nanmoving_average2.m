function Y = nanmoving_average2(X,varargin)
%NANMOVING_AVERAGE2   Smooths a matrix through the moving average method.
%   Y = NANMOVING_AVERAGE(X,F,G) Quickly smooths the matrix X via averaging 
%   each element with their surrounding elements that fit in the little 
%   matrix (2F+1)x(2G+1) centered at the element (boxcar filter), ignoring 
%   the NaN's. The elements at the ends are also averaged but the ones on 
%   the edges are left intact. If F(G) is zero the smoothing is through
%   columns(rows) only. The program uses NANMOVING_AVERAGE. 
%
%   Y = NANMOVING_AVERAGE(X,F), uses G=F.
%
%   Y = NANMOVING_AVERAGE(X), uses G=F=1 by default, it means a 3x3 box.
%
%   Example:
%      [X,Y] = meshgrid(-2:.2:2,3:-.2:-2);
%      Zi = 5*X.*exp(-X.^2-Y.^2); 
%      Zr = Zi + rand(size(Zi)); Zr([8 46 398 400]) = NaN;
%      Zs = nanmoving_average2(Zr,2,3);
%       subplot(131), surf(X,Y,Zi) 
%       view(2), shading interp, xlabel('Z')
%       subplot(132), surf(X,Y,Zr)
%       view(2), shading interp, xlabel('Z + noise + NaN''s')
%       subplot(133), surf(X,Y,Zs)
%       view(2), shading interp, xlabel('Z smoothed')
%
%
%   See also MOVING_AVERAGE, MOVING_AVERAGE2, NANMOVING_AVERAGE

%   Written by
%   M.S. Carlos Adrián Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico, october 2006
%
%   nubeobscura@hotmail.com

% Matrix?
if ndims(X) ~= 2
 disp('ERROR: the entry is not a matrix!')
 return
end
[M,N] = size(X);
Y = zeros(M,N);
suavenans = 0;

% Check entries:
[m,n,suavenans] = checa_arg(varargin,nargin,suavenans);

% Matrix with zeros at NaN elements of X, ones otherwise:
A = double(~isnan(X));

% 1. Sums through columns:
for j = 1:N
 Y(:,j) = nanmoving_average(X(:,j),n,'interpNaN');
 B(:,j) = moving_average(A(:,j),n); % number of elements
end
Y = Y.*B;   % Only sums, do not average

% 2. Smoothing through the sumed rows:
for i = 1:M
 Y(i,:) = nanmoving_average(Y(i,:),m,'interpNaN');
 B(i,:) = moving_average(B(i,:),m);
 C(i,:) = moving_average(A(i,:),m);
end
Y = Y.*C;          % Only sums, do not average
B(B==0) = 1;       % NaN's keep being NaN
Y = Y./B;          % Averaging

% Smooths NaN's?
if ~suavenans
 Y(isnan(X)) = NaN;
end


function [m,n,suavenans] = checa_arg(entradas,Nentr,suavenans)
% Check arguments

if Nentr == 1
 m = 1; n = m; return
end

m = -1; n = -1;
for i = 1:Nentr-1
 if ischar(entradas{i}) && strcmpi(entradas{i}(1),'i')
  suavenans = 1;
 elseif m<0
  m = entradas{i}; 
 elseif n<0
  n = entradas{i};
 end
end
if m<0, m = 1; end
if n<0, n = m; end


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com