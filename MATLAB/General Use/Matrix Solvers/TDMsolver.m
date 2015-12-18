function [x] = TDMsolver(A,b)
%   Solves tridiagonal matrix using Thomas' algorithm.
%Inputs:
%       A:      Full matrix, where Ax=b
%       b:      Column vector
%Outputs:
%       x:      Solution to the linear equation
    
    n = length(b); 
    x = zeros(n,1);
    
    aL = [0; diag(A,-1)];
    a = diag(A);
    aR = [diag(A,1); 0];

    aR(1) = aR(1)/a(1);
    b(1) = b(1)/a(1);
        
    for i=2:n
       t = 1/(a(i)-aR(i-1)*aL(i));
       aR(i) = aR(i)*t;
       b(i) = (b(i)-b(i-1)*aL(i))*t;
    end

    x(n) = b(n);
    for i = n-1:-1:1
       x(i) = b(i)-aR(i)*x(i+1); 
    end    
end