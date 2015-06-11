function [wavelet,fcoef] = MotherWavelets(dim,mother,m)
%This function returns a function handle, 'wavelet', given a dimension and
%string for the name of the mother wavelet. A parameter to modify the base
%mother wavelet can also be provided. For standard mother wavelets, the
%dimension provided will be an integer, for Spatio-Temporal mother
%wavelets, the dimension will be a string beginning with 'ST' and
%including the total number of dimensions.

    if ~exist('m','var'), m = []; end

    switch dim
        case 1
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

        case 'ST2'
            switch lower(mother)
                case 'morlet'
                    if isempty(m), m = [6 6]; end
                    fcoef = m + sqrt(2 + m.^2);
                    wavelet = @(omega,k) (pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k > 0) + ... %positive wavenumbers
                                            (pi^-0.5)*exp(-0.5*((k + m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k < 0); %negative wavenumbers 
                otherwise
                    error('Undefined Mother Wavelet');
            end
        otherwise
            error('Incorrect Dimension Definition');
    end
end