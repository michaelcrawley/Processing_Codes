function [ER,alpha] = WaveletFilterEnergyRatio(sig,dt,mother,param,J,da,amax)
    
    if ~exist('mother','var') || isempty(mother), mother = 'morlet'; end %mother wavelet
    if ~exist('param','var'), param = []; end %wavelet parameter
    if ~exist('J','var') || isempty(J), J = 300; end %number of wavelet scales
    if ~exist('da','var') || isempty(da), da = 0.00005; end %cutoff spacing
    if ~exist('amax','var') || isempty(amax), amax = 0.0015; end %maximum normalized energy cutoff
    
    alpha = 0:da:amax;
    N = size(sig);
    L = 1/(2*dt*N(1));
    ER = zeros(length(alpha),N(2));
    for k = 1:N(3)
        for n = 1:N(2) %step through all channels
            %normalize & remove mean
            x = sig(:,n,k);
            xn = detrend(x)/std(x);

            [wave,~,s,omegak] = CWT1D(xn,dt,1/dt,L,J,mother,param);
            energy = abs(wave*sqrt(dt)).^2; %normalized energy
            for q = 1:length(alpha)
                chk = energy > alpha(q);
                filt = wave;
                filt(~chk) = 0;
                rec = iCWT1D(filt,s,omegak,mother,param);
                ER(q,n) = ER(q,n) + std(rec)*100/N(3); %convert to percent
            end
        end
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