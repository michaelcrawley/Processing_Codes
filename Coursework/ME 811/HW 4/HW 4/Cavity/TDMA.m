function [x] = TDMA(a,b,c,d)
% Tridiagonal matrix algorithm based on the method presented at
% http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
% 	a,b,c are the sub-, main-, and super-diagonals respectively
%	d is the solution vector

N = length(b);
if ~( length(a)==N && length(c)==N && length(d)==N )
	% Check for inproper inputs
	return;
end

% Zero non-existant cells
a(1) = 0; c(N) = 0;

%% Forward sweep
% Normalize the first row
c(1) = c(1)/b(1);
d(1) = d(1)/b(1);

% Sweep through the remaining rows
for n = 2:N
    z = b(n) - c(n-1)*a(n);
    c(n) = c(n)/z;
    d(n) = ( d(n) - d(n-1)*a(n) )/z;
end

clear z a b

%% Backward substitution
% Solve for last row
x(N) = d(N);

% Reverse sweep through the remaining rows
for n = N-1:-1:1
    x(n) = d(n) - c(n)*x(n+1);
end
