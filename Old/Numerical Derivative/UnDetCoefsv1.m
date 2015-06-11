function [coefs] = UnDetCoefsv1(m,n,h,deriv)
    %%Inputs:
    %%m: left side (negative)
    %%n: right side (positive)
    %%h: step size
    %%deriv: derivative to be approximated (1, 2, ...)
    l = n-m+1;
    A = zeros(l,l);
    F = zeros(l,1);
    for i = 0:l-1
        A(:,i+1) = ((m+i)*h).^(0:l-1);
    end
    F(deriv+1) = factorial(deriv);
    coefs = A\F;
end