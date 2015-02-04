function [avg] = wmavg(input,Np,wndo,dim)
    %1 dimensional weighted moving average function
    %avg = wmavg(input,Np,wndo,dim)
    %Inputs:
    %           input: signal to average
    %           Np: number of stensil points (must be odd)
    %           wndo: window for weighting (default is @rectwin)
    %           dim: dimension for averaging (for 2-dimensional arrays)
    %Outputs:
    %           avg: averaged signal

    if ~exist('wndo','var'), wndo = @rectwin; end
    if mod(Np,2) == 0, Np = Np+1; end
    if ~exist('dim','var'), dim = 1; end
    
    weight = window(wndo,Np)';
    
    [M N] = size(input);
    if dim == 2
       input = input.';
       M = N;
    end
    
    x = (Np-1)/2;
    coefs = repmat(weight,M,1);
    for n = 1:x
        scale = sum(coefs(n,1:x+n));
        coefs(n,:) = coefs(n,:)/scale;
        coefs(end+1-n,:) = coefs(end+1-n,:)/scale;
    end
    scale = sum(weight);
    coefs(x+1:end-x,:) = coefs(x+1:end-x,:)/scale;
    
    multiplier = spdiags(coefs,-x:x,M,M)';
    
    avg = multiplier*input;
    if dim == 2
       avg = avg.'; 
    end
end