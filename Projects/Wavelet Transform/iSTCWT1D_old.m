function [xn] = iSTCWT1D(wave,scales,omegaks,mother,param)



    %Integrate over scales - take real part only
    wave = real(wave);
    N = size(wave);
    [a,c] = meshgrid(scales{2},scales{1}); %meshgrid assumes inputs are DIM2,DIM1
    a = repmat(a,[1 1 N(3:4)]);
    c = repmat(c,[1 1 N(3:4)]);
    wave = wave./a./a./c; %normalize by spatial scale and phase velocity
    xn = trapz(scales{2},wave,2); %integrate along spatial scales
    xn = trapz(scales{1},xn,1); %integrate along phase velocities
    xn = squeeze(xn);
    
    %Normalize by C_delta
    wavelet = wavebases(mother,param);
    [K,T] = meshgrid(omegaks{2},omegaks{1}); %order is flipped because meshgrid assumes input is DIM2,DIM1
    Wn = zeros(length(scales{1}),length(scales{2}));
    for n = 1:length(scales{1})
        for m = 1:length(scales{2})
            daughter = wavelet(scales{2}(m)*T/sqrt(scales{1}(n)),scales{2}(m)*K*sqrt(scales{1}(n)));
            Wn(n,m) = mean(real(daughter(:)))/scales{1}(n)/scales{2}(m);
        end
    end
    C_delta = trapz(scales{1},trapz(scales{2},Wn,2)); %Integrate over both dimensions
    
    %Normalize reconstructed signal
    xn = xn/C_delta;
end

function [wavelet,fcoef] = wavebases(mother,m)
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = [6 6]; end
            fcoef = m + sqrt(2 + m.^2);
%             wavelet = @(omega,k) (pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2));%.*(omega < 0);
            wavelet = @(omega,k) (pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k > 0) + ... %positive wavenumbers
                                    (pi^-0.5)*exp(-0.5*((k + m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k < 0); %negative wavenumbers 
        otherwise
            error('Undefined Mother Wavelet');
    end
end