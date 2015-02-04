function [A,recon] = SpectralLSE(field,signal,alpha)
%DIM1 is time
%DIM2 is space
%DIM3 is block
%alpha is the coherence threshold

    if ~exist('alpha','var')||isempty(alpha), alpha = 0; end
    
    tsig = fft(signal,[],1);
    tfield = fft(field,[],1);
    
    [Nfreq,NS,NB] = size(tsig);
    [~,NF,~] = size(tfield);
    
    A = complex(zeros(Nfreq,NS,NF));
    recon = complex(zeros(Nfreq,NF,NB));
    for n = 1:Nfreq %looping through all frequencies
        W = complex(zeros(NS));
        V = complex(zeros(NF,NS));
        svv = complex(zeros(NF,1));
        
        for q = 1:NB
            W = W + tsig(n,:,q)'*tsig(n,:,q);
            V = V + tfield(n,:,q)'*tsig(n,:,q);
            svv = svv + (tfield(n,:,q).*conj(tfield(n,:,q))).';
        end
        W = W/NB;
        V = V/NB;
        svv = svv/NB;
               
        A(n,:,:) = W\V.';   
        
        %Compute coherence and apply to linear weight vectors
        sjj = diag(W);
        [Svv,Sjj] = ndgrid(svv,sjj);
        coherence = (abs(V).^2)./Svv./Sjj;
        chk = coherence >= alpha;
        A(n,:,:) = squeeze(A(n,:,:)).*chk';
    
        for q = 1:NB
            recon(n,:,q) = tsig(n,:,q)*squeeze(A(n,:,:));
        end
    end
    
    recon = ifft(recon,[],1);    
end