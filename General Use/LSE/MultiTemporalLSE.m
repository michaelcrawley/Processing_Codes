function [recon,A] = MultiTemporalLSE(field,signal) 
    
    [BS,NS,NB] = size(signal);    
    NF = size(field,2);
    
    %zero-pad in time, so that correlations can be computed over all time
    %lags
    signal = [zeros(BS-1,NS,NB);signal;zeros(BS-1,NS,NB)];
    field = [zeros(BS-1,NF,NB);field;zeros(BS-1,NF,NB)];
    NC = 2*(3*BS-2)-1;

    conditional = zeros(NC*NS,NF);
    unconditional = zeros(NC*NS);
    for n = 1:NS
        for k = -BS:BS
            
        end
    end

    A = unconditional\conditional;
end