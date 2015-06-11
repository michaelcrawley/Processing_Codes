function L = makeLine(x1,x2,N)
%This function interpolates N points between any two n-dimensional vectors.
%The result is an N by n matrix of points beginning with x1 and ending at
%x2.
%
%INPUTS:
% x1 - starting vector
% x2 - ending vector
% N - the number of points including x1 and x2 to include in the
% interpolation.
%
%OUTPUT:
% L - the matrix of N points arranged as N by n.

x1 = x1(:)';
x2 = x2(:)';

d = x2-x1;

L = repmat(d,N,1).*repmat(linspace(0,1,N)',1,size(x1,2)) +repmat(x1,N,1);