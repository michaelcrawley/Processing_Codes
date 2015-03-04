function [recon,A,f] = SpectralLSE(field,FS_f,signal,FS_s,alpha)
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
%           FS_f:       Sampling frequency (relative or absolute) for field
%           signal:     Generalized reference signal
%           FS_s:       Sampling frequency (relative or absolute) for sig
%           alpha:      Coherence (gamma^2) threshold
%Outputs:
%           recon:      Reconstructed field per signal
%           A:          Frequency-dependent correlation matrix
%           f:          Frequency axis
%Last updated: 2015-03-04 by Michael Crawley

    if ~exist('FS_f','var')||isempty(FS_f), FS_f = 1; end
    if ~exist('FS_s','var')||isempty(FS_s), FS_s = 1; end
    if ~exist('alpha','var')||isempty(alpha), alpha_flag = false; else alpha_flag = true; end
    
    %Get unique frequencies
    [nfft_s,NS,NB] = size(signal);
    [nfft_f,NF,~] = size(field);
    nfrq_s = ceil((nfft_s+1)/2);
    nfrq_f = ceil((nfft_f+1)/2);
    f_s = (0:nfrq_s-1)*(FS_s/nfft_s); 
    f_f = (0:nfrq_f-1)*(FS_f/nfft_f);
    f = union(f_s,f_f); %full frequency axis - includes both common and non-common frequencies
    [uf,ifs,iff] = intersect(f_s,f_f); %common frequencies
    [~,iuf] = intersect(f,uf); %common frequencies (indices) only
    nfrq = length(f);
    
    %compute fft
    if mod(nfrq,2) %nyquist frequency is measured - nonsymmetric positive/negative frequency vector
        nfft = (nfrq-1)*2;
    else
        nfft = (nfrq-1)*2+1;
    end
    tsig = fft(signal,[],1)*(nfft/nfft_s); %fix energy scaling
    tfield = fft(field,[],1)*(nfft/nfft_f); %fix energy scaling
    
    A = complex(zeros(NS,NF,nfrq));
    recon = complex(zeros(nfrq,NF,NB));
    for n = 2:length(iuf) %looping through all relevant frequencies - skip zero frequency
        
        W = complex(zeros(NS)); %unconditional signal
        V = complex(zeros(NF,NS)); %conditional field
        if alpha_flag
            svv = complex(zeros(NF,1));
        end
        
        for q = 1:NB
            W = W + tsig(ifs(n),:,q)'*tsig(ifs(n),:,q);
            V = V + tfield(iff(n),:,q)'*tsig(iff(n),:,q);
            if alpha_flag
                svv = svv + (tfield(iff(n),:,q).*conj(tfield(iff(n),:,q))).';
            end
        end
        W = W/NB;
        V = V/NB;
        if alpha_flag
            svv = svv/NB;
        end
               
        chk = W == 0;
        if all(chk(:)), W = Inf; end
        A(:,:,iuf(n)) = W\V';   
        
        %Compute coherence and apply to linear weight vectors
        if alpha_flag
            sjj = diag(W);
            [Svv,Sjj] = ndgrid(svv,sjj);
            coherence = (abs(V).^2)./Svv./Sjj;
            chk = coherence >= alpha;
            A(:,:,iuf(n)) = A(:,:,iuf(n)).*chk;
        end
    
        for q = 1:NB
            recon(iuf(n),:,q) = tsig(ifs(n),:,q)*A(:,:,iuf(n));
        end
    end
    A = permute(A,[3 1 2]); %switch so that frequency is first dimension - its easier to concatenate this way
    
    %For the iFFT, both the positive and negative frequencies are necessary
    if mod(nfft,2) %signal length is odd
        A = [A;conj(flipud(A(2:end,:,:)))];
        recon = [recon;conj(flipud(recon(2:end,:,:)))];
    else %signal length is even
        inyquist = min(nfrq_f,nfrq_s);
        A(inyquist,:,:) = A(inyquist,:,:)/2;
        recon(inyquist,:,:) = recon(inyquist,:,:)/2; %halve energy at nyquist frequency
        A = [A;conj(flipud(A(2:end-1,:,:)))];
        recon = [recon;conj(flipud(recon(2:end-1,:,:)))];
    end
    
    recon(1,:,:) = tfield(1,:,:); %insert mean value
    recon = ifft(recon,[],1);    
end