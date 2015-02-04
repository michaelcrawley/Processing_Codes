function [coefs] = nuTSE(sp,x,deriv)
    %Calculate matrix coefficients for use in finite difference by Taylor series expansions.
    %Code Version: 3.0 @ 2015-01-28
    %Inputs:
    %   sp:     number of stencil points (should be odd number)
    %   x:      grid points 
    %   deriv:  derivative to be approximated (1, 2, ...)
    %Output:
    %   coefs:  
    
    if ~mod(sp,2)
        warning('Switching to odd-numbered stencil for central differencing...');
        sp = sp+1;
    end
    coefs = zeros(length(x));

    F = zeros(sp,1);
    F(deriv+1) = factorial(deriv);
    stencil = -(sp-1)/2:(sp-1)/2;
    
    for n = 1:length(x)
        loc = stencil + n;
        if any(loc <= 0), loc = loc + abs(loc(1)) + 1; end
        if any(loc > length(x)), loc = loc - (max(loc) - length(x)); end
        h = x(loc) - x(n);
        A = zeros(sp);
        for q = 1:sp
            A(q,:) = h.^(q-1);
        end
        coefs(n,loc) = (A\F);
    end    
end