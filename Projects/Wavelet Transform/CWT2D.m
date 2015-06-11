function [wave,freq,s,omegak] = CWT2D(sig,dt,dj,mother,param)
%Computes 2-Dimensional Continuous Wavelet Transform using Fourier method.
%Inputs:
%           x:          signal (M,N)
%           dt:         period between samples (2,1)
%           dq:         scale spacings (2,1)
%           mother:     mother wavelet (Default is Morlet)
%           param:      parameter for mother wavelet
%Outputs:
%           wave:       wavelet coefficients (J,Q,M,N)
%           freq:       Fourier-equivalent frequency vector
%           s:          scale index used {2,1}
%           omegak:     frequency index used for FFT {2,1}
%Last Updated: 2014-04-14 by Michael Crawley
    

    %Grab Mother Wavelet info
    if ~exist('mother','var'), mother = 'morlet'; end
    if ~exist('param','var'), param = []; end
    if ~exist('dj','var'), dj = [0.25 0.25]; end
    [wavelet,fcoef] = wavebases(mother,param);  
    
    %Get spacings
    [M,N] = size(sig);
    s0 = 2*dt(1)*fcoef(1);
    c0 = 2*dt(2)*fcoef(2);
    J = round(log2(M*dt(1)/s0)/dj(1))+1;
    Q = round(log2(N*dt(2)/c0)/dj(2))+1;
    s{1} = s0*2.^(dj(1)*(0:J-1)); %DIM1 scale
    s{2} = c0*2.^(dj(2)*(0:Q-1)); %DIM2 scale

    %Compute FFT Processing parameters
    omega = 2*pi*(1:fix(M/2))/M/dt(1);
    omegak{1} = [0,omega,-fliplr(omega(1:fix((M-1)/2)))]; %FFT frequency index - DIM1
    omega = 2*pi*(1:fix(N/2))/N/dt(2);
    omegak{2} = [0,omega,-fliplr(omega(1:fix((N-1)/2)))]; %FFT frequency index - DIM2
    [Wk2,Wk1] = meshgrid(omegak{2},omegak{1});
    sigT = fftn(sig);
    freq = {1./(fcoef(1)*s{1}),1./(fcoef(2)*s{2})};
    
    %Initialize output, loop through scale index
    wave = zeros(J,Q,M,N);
    for j = 1:J
        for q = 1:Q
            daughter = wavelet(s{1}(j)*Wk1, s{2}(q)*Wk2);
            wave(j,q,:,:) = ifftn(sigT.*conj(daughter)*((2*pi*s{1}(j)/dt(1))^0.5)*((2*pi*s{2}(q)/dt(2))^0.5));
        end
    end
end

function [wavelet,fcoef] = wavebases(mother,m)
    switch lower(mother)
        case 'morlet'
            if isempty(m), m = [6 6]; end
            wavelet = @(omega,k) (pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(omega > 0); 
            fcoef = (4*pi)./(m + sqrt(2 + m.^2)); %k_m and omega_m 
        otherwise
            error('Undefined Mother Wavelet');
    end
end