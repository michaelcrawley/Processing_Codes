function [coh phase cpsd xpsd ypsd] = WaveletCoherence(x,y,scale,dt)
    %Computes the wavelet coherence between wavelet coefficients for two
    %signals. Currently, this function only works for wavelet coefficients
    %computed using the morlet mother wavelet, as the time
    %smoothing function for other wavelets are not currently known by me.
    %This function is based off of Torrence 1998 & 1999.
    %Inputs:
    %           x:      Wavelet Coefficients (complex) for signal 1 (MxN matrix)
    %                   where M is number of scales and N is number of time steps 
    %           y:      Wavelet Coefficients for signal 2 (MxN matrix)
    %           scale:  Vector of scales used (Mx1 matrix)
    %           dj:     constant used when computing scale vector
    %           dt:     sampling time of signals
    %Outputs:
    %           coh:    Wavelet-Coherence map (MxN matrix)
    %           phase:  Wavelet-Coherence phase map (MxN matrix)
    
    %Create Smoothing Functions
    [M,N] = size(x); 
    t = ((1:N)-floor(N/2))*dt; %time series, centered
    time_smooth = @(s) exp(-(t.^2)/(2*s^2)).'; %smoothing window for morlet wavelet only!
    dj = 0.6; %scale decorrelation length for morlet wavelet only!
    for n = 1:M
        ssm(n,:) = abs(scale(n)-scale) < dj/2*scale(n);
    end
    ssm = ssm./repmat(sum(ssm,2),1,M); %normalize by convolution energy at each scale -- Is this correct???
    
    %Compute wavelet power spectra and cross spectra - for simplicity,
    %these will not be properly normalized (as coherence is by definition
    %normalized)
    xpsd = abs(x).^2; %absolute values for PSD's are taken before smoothing, after smoothing for cross-PSD
    ypsd = abs(y).^2;
    cpsd = x.*conj(y);
    
    %Smooth in time at each scale
    for n = 1:M
        tsv = time_smooth(scale(n));
        tsm = spdiags(repmat(tsv,1,N),(1:N)-floor(N/2)); %smoothing matrix
        tsm = tsm./repmat(sum(tsm,1),N,1); %normalize by convolution energy at each time step -- Is this correct????
        
        %smooth
        xpsd(n,:) = xpsd(n,:)*tsm;
        ypsd(n,:) = ypsd(n,:)*tsm;
        cpsd(n,:) = cpsd(n,:)*tsm;
    end
    
    %Smooth in scale at each time step
    xpsd = ssm*xpsd;
    ypsd = ssm*ypsd;
    cpsd = ssm*cpsd;
    
    %Compute Coherence and Phase
    coh = abs(cpsd).^2./(xpsd.*ypsd); %need to fix energy in cross-spectra 
    phase = atan(imag(cpsd)./real(cpsd));
end