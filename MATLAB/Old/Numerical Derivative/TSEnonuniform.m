function [coefs] = TSE(n,m,x,deriv)
    %Calculate coefficients for use in finite difference by Taylor series expansions.
    %Code Version: 2.1 @ 2011-04-03
    %Inputs:
    %   n: number of backwards terms 
    %   m: number of forwards terms (stencil goes from -n:m)
    %   x: either step size (scalar) or grid points (array)
    %   deriv: derivative to be approximated (1, 2, ...)

    l = n+m+1;
    A = zeros(l,l);
    F = zeros(l,1);
    if length(x) == 1
       h = x*(-n:m);
    else
       h = x-x(n+1);
    end
    
    for i = 1:l
        A(i,:) = h.^(i-1);
    end

    F(deriv+1) = factorial(deriv);
    coefs = (A\F);
end