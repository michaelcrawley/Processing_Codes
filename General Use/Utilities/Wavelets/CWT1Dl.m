function [wave,freq,s] = CWT1Dl(sig,dt,fmax,fmin,J,mother,param)
%Computes 1-Dimensional Continuous Wavelet Transform using a standard for-loop Convolution method.
%Inputs:
%           x:          signal (N,1)
%           dt:         period between samples
%           fmax:       maximum Fourier frequency of interest
%           fmin:       minimum Fourier frequency of interest
%           J:          number of scales to use
%           mother:     mother wavelet (Default is Morlet)
%           param:      parameter for mother wavelet
%Outputs:
%           wave:       wavelet coefficients (J,N)
%           freq:       Fourier-equivalent frequency vector
%           s:          scale index used
%Last Updated: 2014-09-02 by Michael Crawley
    

    %Grab Mother Wavelet info
    if ~exist('mother','var'), mother = 'morlet'; end
    if ~exist('param','var'), param = []; end
    [wavelet,fcoef] = wavebases(mother,param);  
    so = 1/(fcoef*fmax);
    sL = 1/(fcoef*fmin);

    %Compute Processing parameters
    N = length(sig);
    t = (1:N)*dt;
    dj = log2(sL/so)/(J-1);
    s = so*2.^(dj*(0:J-1)); %scale index
    freq = 1./(fcoef*s);
    
    %Initialize output, loop through scale index
    wave = zeros(J,N);
    for j = 1:J
        for n = 1:N
            daughter = wavelet((n*dt-t)/s(j))*(dt/s(j))^0.5;
%             daughter = daughter/sqrt(trapz(abs(daughter).^2)); %norm to remove edge effects
            wave(j,n) = trapz(sig.*conj(daughter));
        end
    end
end

function [wavelet,fcoef] = wavebases(mother,m)
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = 6; end
            fcoef = (4*pi)/(m + sqrt(2 + m^2));
            wavelet = @(x) (pi^-0.25)*exp(1i*m*x).*exp(-(x.^2)/2); 
        otherwise
            error('Undefined Mother Wavelet');
    end
end