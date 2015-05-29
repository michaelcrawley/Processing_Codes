function [ N, M, S ] = nzstats( x, dim )
%NZSTATS Calculate the mean and standard deviation of nonzero elements.

% Number of non-zero vectors
N = sum( x~=0, dim );

D = max( N, 1 );			% Avoid division by zero

% Mean
M = sum( x, dim ) ./ D;

% Standard deviation
xf = bsxfun( @minus, x, M );
S = sqrt( sum( xf.^2, dim ) ./ D );