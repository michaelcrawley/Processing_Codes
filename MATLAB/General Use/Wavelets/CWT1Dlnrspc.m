function [wave,freq,s,omegak] = CWT1Dlnrspc(x,dt,fmax,fmin,J,mother,param)
%Computes 1-Dimensional Continuous Wavelet Transform using Fourier method.
%For this code, the scales are linearly spaced, specifically for use with
%POD (so I don't need to apply a weighting function to the grid).
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
%           omegak:     Normalized Fourier frequencies for FFT of x
%Last Updated: 2014-03-18 by Michael Crawley
    

    %Grab Mother Wavelet info
    if ~exist('mother','var'), mother = 'morlet'; end
    if ~exist('param','var'), param = []; end
    [wavelet,fcoef] = wavebases(mother,param);  
    so = 1/(fcoef*fmax);
    sL = 1/(fcoef*fmin);

    %Compute Processing parameters
    N = length(x);
    omega = 2*pi*(1:floor(N/2))/N/dt;
    omegak = [0,omega,-fliplr(omega(1:floor((N-1)/2)))]; %FFT frequency index
%     dj = log2(sL/so)/(J-1);
%     s = so*2.^(dj*(0:J-1)); %scale index
    s = linspace(so,sL,J);
    xh = fft(x);
    freq = 1./(fcoef*s);
    
    %Initialize output, loop through scale index
    wave = zeros(J,N);
    for j = 1:J
        daughter = wavelet(s(j)*omegak).';
        wave(j,:) = ifft(xh.*conj(daughter))*(2*pi*s(j)/dt)^0.5;        
    end
end

function [wavelet,fcoef] = wavebases(mother,m)
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = 6; end
            fcoef = (4*pi)/(m + sqrt(2 + m^2));
            wavelet = @(k) (pi^-0.25).*exp(-0.5*((k - m).^2)).*(k > 0); 
        case 'paul'
            if isempty(m), m = 4; end
            fcoef = 4*pi/(2*m+1);
            wavelet = @(k) (2^m/sqrt(m*factorial(2*m-1))).*(k.^m).*exp(-(k).*(k > 0)).*(k > 0);
        case 'dog'
            if isempty(m), m = 2; end
            fcoef = 2*pi/sqrt(m+0.5);
            wavelet = @(k) -(1i^m)/sqrt(gamma(m+0.5))*(k.^m).*exp(-0.5*k.^2);
        otherwise
            error('Undefined Mother Wavelet');
    end
end
