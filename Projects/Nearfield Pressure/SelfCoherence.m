function [coh f xpsd] = SelfCoherence(x,FS,BS,whan)
%Computes the coherence and phase between two signals
%[f coh phase cpsd xpsd ypsd] = coherence(x,y,FS,BS,whan)
%Inputs: 
%           x: first signal (M x N)
%           FS: sampling rate (Hz)
%           BS: block size for computations (does not have to equal M)
%           whan: window handle
%Outputs:
%           coh: estimated self-coherence (M x M)
%           f: frequency axis (Hz)
%           xpsd: power spectral density for x signal
%Last modified on 2014-08-18 by Michael Crawley
    

    %determine array sizes
    [N NB] = size(x);
    if isempty(BS), BS = NB; end
    if N ~= BS
        x = reshape(x,BS,[]);
        [~, NB] = size(x);
    end
    
    %Calculate frequency axis
    nfrq = ceil((BS+1)/2);  %number of distinct frequencies returned for SPL
    f = (0:nfrq-1)'*(FS/BS);   %Frequency axis
       
    %Create window
    wndo = whan(BS);
    ww = mean(wndo.^2); %window weigth
    wndo = repmat(wndo,1,NB);
    
    %normalize signal
    xm = mean(x,1);
    xn = x - repmat(xm,BS,1);
    
    %compute FFT
    xt = fft(xn.*wndo);
    
    %get rid of symmetric parts
    xt = xt(1:nfrq,:);
    xt(2:end-1) = xt(2:end-1)*sqrt(2);
    
    %compute power spectral densities
    xpsd = mean(conj(xt).*xt,2)/(BS*FS*ww);
     
    %compute coherence
    cpsd = (xt*xt')/(BS*FS*ww)/NB; %note that I am intentionally taking the conjugate transpose
    coh = sqrt(abs(cpsd).^2./(xpsd*xpsd.'));
end