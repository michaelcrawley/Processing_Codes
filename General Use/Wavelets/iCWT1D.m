function [xn] = iCWT1D(wave,scale,omegak,mother,param)
%Computes inverse 1-Dimensional Continuous Wavelet Transform using Fourier 
%method. This algorithm performs the integration over scale, rather than
%over translation (due to cone-of-influence effects near the boundaries).
%Inputs:
%           wave:       wavelet coefficients (M,N)
%           scale:      scale vectors (M,1)
%           omegak:     Normalized Fourier frequencies for FFT of x
%           mother:     mother wavelet 
%           param:      parameter for mother wavelet
%Outputs:
%           xn:         Reconstructed signal
%Last Updated: 2014-04-24 by Michael Crawley

    [M,N] = size(wave);
    if M ~= length(scale), wave = wave.'; N = M; end
    wavelet = wavebases(mother,param);
    
    %Integrate wavelet coefficients using trapezoidal method
    s = repmat(scale(:),[1,N]).^1.5;
    real_norm_waves = real(wave)./s;
    xn = trapz(scale,real_norm_waves,1); %integrate along scale dimension
    
    %Calculate C_delta for normalization
    Wn = zeros(length(scale),1);
    for n = 1:length(scale)
        Wn(n) = sqrt(omegak(2)*N)*mean(real(wavelet(scale(n)*omegak)))/(scale(n)^1); 
    end
    C_delta = trapz(scale,Wn);
    
    %Normalize reconstructed signal
    xn = xn/C_delta;
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