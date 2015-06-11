function [coh f cpsd xpsd ypsd] = coherence_v2(x,y,FS,BS,whan)
%Computes the coherence and phase between two signals. This code is
%designed to check the results of coherence.m. In this code, the magnitude
%of the cross-spectral density is computed prior to the block-averaging, as
%opposed to afterwards as it is done in coherence.m.
%[f coh phase cpsd xpsd ypsd] = coherence(x,y,FS,BS,whan)
%Inputs: 
%           x: first signal (M x N)
%           y: second signal (M x N)
%           FS: sampling rate (Hz)
%           BS: block size for computations (does not have to equal M)
%           whan: window handle
%Outputs:
%           f: frequency axis (Hz)
%           coh: estimated coherence (gamma ^ 2)
%           cpsd: cross-power spectral density
%           xpsd: power spectral density for x signal
%           ypsd: power spectral density for y signal
    

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
    ypsd = 2*mean(conj(yt).*yt,2)/(BS*FS*ww);
    cpsd = 2*mean(abs(conj(xt).*yt),2)/(BS*FS*ww);
     
    %compute coherence (gamma^2) and phase
    coh = cpsd.^2./(xpsd.*ypsd);
end
