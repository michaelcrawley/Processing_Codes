function [x] = iPSDN(S,order,xm)
   
    %Compute FFT
    N = length(order);
    x = S;
    for n = 1:N
        x = ifft(x,[],order(n));
    end
    
    %Reinsert mean
    x = x+xm;
end