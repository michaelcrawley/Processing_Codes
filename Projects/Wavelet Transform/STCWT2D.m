function [wave,s,c] = STCWT2D(sig,dt,dx,dj,dq,mother,param)
%Computes 1-Dimensional Continuous Wavelet Transform using Fourier method.
%%%%%%%%%%%%%%%%[UNFINISHED]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%           sig:          signal (M,N)
%           dt:         period between M-samples
%           dx:         period between N-samples
%           dj:         scale spacing
%           dq:         phase velocity spacing
%           mother:     mother wavelet (Default is Morlet)
%           param:      parameter for mother wavelet
%Outputs:
%           wave:       wavelet coefficients (J,Q,M,N)
%           s:          scale index used
%           c:          phase velocity index used
%Last Updated: 2014-04-11 by Michael Crawley    

    %Grab Mother Wavelet info
    if ~exist('mother','var'), mother = 'morlet'; end
    if ~exist('param','var'), param = []; end
    [wavelet,coefs] = wavebases(mother,param);  
    
    %Get spacings
    [M,N] = size(sig);
    s0 = 1/(2*pi)*sqrt(prod(coefs))*sqrt(dt*dx);
    J = round((log2(sqrt(M*N))-1)/dj)+1;
    s = s0*2.^(dj*(0:J-1));
    c0 = 2*dx*coefs(1)/(M*dt*coefs(2));
    Q = round((log2(M*N)-2)/dq)+1;
    c = c0*2.^(dq*(0:Q-1));
    
end

function [wavelet,coefs] = wavebases(mother,m)
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = [6 6]; end
            wavelet = @(k,omega) (pi^-0.125)*exp(-0.5*((k - m(1)).^2)).*exp(-0.5*((omega - m(2)).^2)).*(omega > 0); 
            coefs = m+sqrt(m.^2+2); %k_m and omega_m
        otherwise
            error('Undefined Mother Wavelet');
    end
end