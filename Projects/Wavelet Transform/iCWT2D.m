function [xn] = iCWT2D(wave,scale,omegak,mother,param)

    %Grab mother wavelet
    mother = wavebases(mother,param);

    %Take real part of wavelet coefficients, normalize by scales, and sum
    N = size(wave);
    wave = real(wave);
    [S2,S1] = meshgrid(sqrt(scale{2}),sqrt(scale{1}));
    S1 = repmat(S1,[1 1, N(3:4)]);
    S2 = repmat(S2,[1 1, N(3:4)]);
    xn = squeeze(sum(sum(wave./S1./S2,1),2));
    
    %Compute C_delta
    W_delta = zeros(N(1:2));
    [Wk2,Wk1] = meshgrid(omegak{2},omegak{1});
    for n = 1:N(1)
        for m = 1:N(2)
            daughter = mother(scale{1}(n)*Wk1,scale{2}(m)*Wk2);
            W_delta(n,m) = sum(conj(daughter(:)))/prod(N(1:2));
        end
    end
    W_delta = W_delta./S1(:,:,1,1)./S2(:,:,1,1);
    C_delta = sum(real(W_delta(:)));
    
    %Normalize
    xn = xn/C_delta;
end

function [wavelet,fcoef] = wavebases(mother,m)
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = [6 6]; end
            wavelet = @(k,omega) (pi^-0.5)*exp(-0.5*((k - m(1)).^2)).*exp(-0.5*((omega - m(2)).^2)).*(omega > 0); 
            fcoef = (4*pi)./(m + sqrt(2 + m.^2)); %k_m and omega_m 
        otherwise
            error('Undefined Mother Wavelet');
    end
end