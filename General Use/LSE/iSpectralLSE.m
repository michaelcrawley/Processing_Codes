function [recon] = iSpectralLSE(f_c,A,signal,fs)
%The purpose of this function is to perform a conditional reconstruction of
%a generalized field from a generalized signal, per the theory of linear
%stochastic estimation. The frequency-dependent correlation matrix will be
%supplied by SpectralLSE.m. The specific function of this routine is to
%effectively upsample the generalized field using the higher temporal
%resolution of the generalized signal.
%Inputs:
%           f_c:        Frequency axis for correlation matrix
%           A:          Frequency-dependent correlation matrix
%           signal:     Generalized reference signal
%           fs:         Sampling frequency for reference signal
%Outputs:
%           recon:      Reconstructed field
%Last updated: 2015-02-24 by Michael Crawley

    %Calc normalized frequency axes
    nf_c = f_c(:)/f_c(2); %normalize by lowest frequency
    nfft = size(signal,1);
    nfrq = ceil((nfft+1)/2);
    f_s = (0:nfrq-1)*(fs/nfft); %frequency axis for the reference signal
    nf_s = f_s(:)/f_c(2); %normalize by lowest frequency in correlation matrix
        
    %build upsampled correlation matrix
    N = size(A);
    A_up = zeros([nfrq,N(2:end)]);
    for n = 1:length(f_c)
        indx = (1:length(f_s))*(nf_c(n) == nf_s); 
        A_up(indx,:,:) = A(n,:,:);
    end
    
    %perform reconstruction
    tsig = fft(signal,[],1);
    NB = size(signal,3);
%     recon = complex(zeros(Ns(1),N(2),Ns(3)));
    for n = 1:nfrq
        for q = 1:NB %fix later for multiple blocks
            recon(n,:,q) = tsig(n,:,q)*squeeze(A_up(n,:,:));
        end
    end
    
    %create negative frequencies
    if mod(nfft,2) %signal length is odd
        recon = [recon;conj(flipud(recon(2:end,:,:)))];
    else %signal length is even
        recon = [recon;conj(flipud(recon(2:end-1,:,:)))];
    end
    
    recon = ifft(recon,[],1);  
end