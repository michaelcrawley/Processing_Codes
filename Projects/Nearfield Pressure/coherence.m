function [coh f phase cpsd xpsd ypsd] = coherence(x,y,FS,BS,whan)
%Computes the coherence and phase between two signals
%[f coh phase cpsd xpsd ypsd] = coherence(x,y,FS,BS,whan)
%Inputs: 
%           x: first signal (M x N)
%           y: second signal (M x N)
%           FS: sampling rate (Hz)
%           BS: block size for computations (does not have to equal M)
%           whan: window handle
%Outputs:
%           coh: estimated coherence (gamma)
%           f: frequency axis (Hz)
%           phase: estimated phase
%           cpsd: cross-power spectral density
%           xpsd: power spectral density for x signal
%           ypsd: power spectral density for y signal
%Last modified on 2012-11-07 by Michael Crawley
    

    %determine array sizes
    [N NB] = size(x);
    if isempty(BS), BS = NB; end
    if N ~= BS
        x = reshape(x,BS,[]);
        y = reshape(y,BS,[]);
        [~, NB] = size(x);
    end
    
    %Calculate frequency axis
    nfrq = ceil((BS+1)/2);  %number of distinct frequencies returned for SPL
    f = (0:nfrq-1)'*(FS/BS);   %Frequency axis
       
    %Create window
    wndo = whan(BS);
    ww = mean(wndo.^2); %window weigth
    wndo = repmat(wndo,1,NB);
    
    %normalize signals
    xm = mean(x,1);
    ym = mean(y,1);
    xn = x - repmat(xm,BS,1);
    yn = y - repmat(ym,BS,1);
    
    %compute FFT
    xt = fft(xn.*wndo);
    yt = fft(yn.*wndo);
    
    %get rid of symmetric parts
    xt = xt(1:nfrq,:);
    yt = yt(1:nfrq,:);
    
    %compute cross and power spectral densities
    xpsd = 2*mean(conj(xt).*xt,2)/(BS*FS*ww);
    xpsd([1,end]) = xpsd([1,end])/2; %energy at f = 0 and f = Nyquist should not be doubled
    ypsd = 2*mean(conj(yt).*yt,2)/(BS*FS*ww);
    ypsd([1,end]) = ypsd([1,end])/2;
    cpsd = 2*mean(conj(xt).*yt,2)/(BS*FS*ww);
     
    %compute coherence (gamma^2) and phase
    coh = sqrt(abs(cpsd).^2./(xpsd.*ypsd));
    phase = atan2(imag(cpsd),real(cpsd));
end
