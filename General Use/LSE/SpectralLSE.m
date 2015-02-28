function [recon,A,f,c] = SpectralLSE(field,signal,fs,alpha)
%This function computes a least-squares correlation matrix between a given
%field and a reference signal, per the theory outlined by Adrian regarding
%conditional averaging based on linear stochastic estimation. For
%simplicity, the correlations are computed in the Fourier domain, so that
%propagation delays can be neglected. A coherence threshold is included as
%an option in order to remove the influence of low-correlated frequencies.
%The field and reference signals need to be ordered as either 3-D or 2-D
%matrices, the the first dimension being time, the second being
%space/component/etc, and the third being separate blocks (block-averaging
%of the correlation matrix is performed by default). 
%Inputs:
%           field:      Generalized field for conditional averaging
%           signal:     Generalized reference signal
%           fs:         Sampling frequency (optional)
%           alpha:      Coherence (gamma^2) threshold
%Outputs:
%           recon:      Reconstructed field per signal
%           A:          Frequency-dependent correlation matrix
%           f:          Frequency axis
%Last updated: 2015-02-24 by Michael Crawley

    if ~exist('fs','var')||isempty(fs), fs = 1; end
    if ~exist('alpha','var')||isempty(alpha), alpha_flag = false; else alpha_flag = true; end
    
    tsig = fft(signal,[],1);
    tfield = fft(field,[],1);
    
    [nfft,NS,NB] = size(tsig);
    [~,NF,~] = size(tfield);
    
    %Get unique frequencies
    nfrq = ceil((nfft+1)/2);
    f = (0:nfrq-1)*(fs/nfft); %frequency axis
    
    A = complex(zeros(NS,NF,nfrq));
    recon = complex(zeros(nfrq,NF,NB));
    for n = 2:nfrq %looping through all relevant frequencies - skip zero frequency
        W = complex(zeros(NS));
        V = complex(zeros(NF,NS));
        if alpha_flag
            svv = complex(zeros(NF,1));
        end
        
        for q = 1:NB
            W = W + tsig(n,:,q)'*tsig(n,:,q);
            V = V + tfield(n,:,q)'*tsig(n,:,q);
            if alpha_flag
                svv = svv + (tfield(n,:,q).*conj(tfield(n,:,q))).';
            end
        end
        W = W/NB;
        V = V/NB;
        if alpha_flag
            svv = svv/NB;
        end
               
        A(:,:,n) = W\V';   
        
        %Compute coherence and apply to linear weight vectors
        if alpha_flag
            sjj = diag(W);
            [Svv,Sjj] = ndgrid(svv,sjj);
            coherence = (abs(V).^2)./Svv./Sjj;
            chk = coherence >= alpha;
            A(:,:,n) = A(:,:,n).*chk;
        end
    
        for q = 1:NB
            recon(n,:,q) = tsig(n,:,q)*A(:,:,n);
        end
    end
    A = permute(A,[3 1 2]); %switch so that frequency is first dimension
    
    %For the iFFT, both the positive and negative frequencies are necessary
    if mod(nfft,2) %signal length is odd
        A = [A;conj(flipud(A(2:end,:,:)))];
        recon = [recon;conj(flipud(recon(2:end,:,:)))];
    else %signal length is even
        A = [A;conj(flipud(A(2:end-1,:,:)))];
        recon = [recon;conj(flipud(recon(2:end-1,:,:)))];
    end
    c = recon;
    
    recon = ifft(recon,[],1);    
end