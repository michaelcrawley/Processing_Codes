function Y = moving_average2(X,varargin)
%MOVING_AVERAGE2   Smooths a matrix through the moving average method.
%   Y = MOVING_AVERAGE(X,F,G) Quickly smooths the matrix X via averaging 
%   each element with their surrounding elements that fit in the little 
%   matrix (2F+1)x(2G+1) centered at the element (boxcar filter). The 
%   elements at the ends are also averaged but the ones on the edges are 
%   left intact. If F(G) is zero the smoothing is over each column(row) 
%   only. The program uses MOVING_AVERAGE. 
%
%   Y = MOVING_AVERAGE(X,F), uses G=F.
%
%   Y = MOVING_AVERAGE(X), uses G=F=1 by default, it means a 3x3 box.
%
%   Example:
%      [X,Y] = meshgrid(-2:.2:2,3:-.2:-2);
%      Zi = 5*X.*exp(-X.^2-Y.^2); 
%      Zr = Zi + rand(size(Zi));
%      Zs = moving_average2(Zr,2,3);
%       subplot(131), surf(X,Y,Zi) 
%       view(2), shading interp, xlabel('Z')
%       subplot(132), surf(X,Y,Zr)
%       view(2), shading interp, xlabel('Z + noise')
%       subplot(133), surf(X,Y,Zs)
%       view(2), shading interp, xlabel('Z smoothed')
%
%
%   See also MOVING_AVERAGE, NANMOVING_AVERAGE, NANMOVING_AVERAGE2

%   Written by
%   M.S. Carlos Adrián Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico, september 2006
%
%   nubeobscura@hotmail.com

% October 2006, fixed bug on the rows.

% Matrix?
if ndims(X) ~= 2
 disp('ERROR: the entry is not a matrix!')
 return
end
[M,N] = size(X);
Y = zeros(M,N);

% Check entries:
if nargin == 1
 m = 1; n = m;
elseif nargin == 2
 m = varargin{1}; n = m;
else
 m = varargin{1}; n = varargin{2};
end

% 1. Smooths through columns:
for j = 1:N
 Y(:,j) = moving_average(X(:,j),n);
end

% 2. Smooths through smoothed rows:
for i = 1:M
 Y(i,:) = moving_average(Y(i,:),m); % Fixed bug
end


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com