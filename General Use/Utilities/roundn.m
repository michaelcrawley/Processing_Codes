function out = roundn(x,n)
%Rounds to the n-th digit
    out = (10^n)*round(x/10^n);
end