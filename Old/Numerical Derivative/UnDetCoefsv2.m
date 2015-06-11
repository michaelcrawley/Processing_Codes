function [coefs] = UnDetCoefsv2(x,h,direc,deriv)
    %Calculate coefficients for use in finite difference by Method of
    %Undetermined Coefficients.
    %Code Version: 1.1 @ 2010-11-14
    %Inputs:
    %   x: number of terms 
    %   h: step size
    %   direc: "center", "down", or "up"
    %   deriv: derivative to be approximated (1, 2, ...)
    if strcmp(direc,'center')
        n = x;
        m = -x;
    elseif strcmp(direc,'down')
        n = x;
        m = 0;
    elseif strcmp(direc,'up')
        n = 0;
        m = -x;
    end
    l = n-m+1;
    A = zeros(l,l);
    F = zeros(l,1);
    for i = 0:l-1
        A(:,i+1) = ((m+i)*h).^(0:l-1);
    end
    F(deriv+1) = factorial(deriv);
    coefs = (A\F);
end