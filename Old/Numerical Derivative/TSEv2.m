function [coefs] = TSE(n,m,h,deriv)
    %Calculate coefficients for use in finite difference by Taylor series expansions.
    %Code Version: 2.0 @ 2011-04-03
    %Inputs:
    %   n: number of backwards terms 
    %   m: number of forwards terms (stencil goes from -n:m)
    %   h: step size (either scalar or array)
    %   deriv: derivative to be approximated (1, 2, ...)

    l = n+m+1;
    A = zeros(l,l);
    F = zeros(l,1);
    if length(h) == 1
       h = h*ones(l,1); 
    end
    
    for i = 0:l-1
        A(:,i+1) = ((-n+i)*h(i+1)).^(0:l-1);
    end

    F(deriv+1) = factorial(deriv);
    coefs = (A\F);
end