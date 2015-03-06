function [recon,A,f] = MultiSpectralLSE(field,signal,FS)

%This is a highly-specialized code, designed specifically to upsample
%non-time resolved PIV data (conditional field) using a generic set of
%microphone data (unconditional signal). This code utilizes the algorithm
%laid out by Tinney 2008 JFM.
%Inputs:
%

    if ~exist('FS','var')||isempty(FS), FS = 1; end
    
    [BS,NF,NB] = size(field);
    NS = size(signal,2);



    %Transform into Fourier domain and compute estimation coefficients at
    %each frequency
    unconditional = fft(unconditional,[],3);
    conditional = fft(conditional,[],3);
    nfrq = 2*Nc+1;
    f = (-Nc:Nc)*(FS/nfrq);
    A = complex(zeros(NS,NF,nfrq));
    for n = 1:nfrq
        A(:,:,n) = unconditional(:,:,n)\conditional(:,:,n);
    end
    A(isnan(A)|isinf(A)) = 0;%fix for badly conditioned sets
    
    %Compute reconstruction in Fourier domain, transform into time domain
    tsig = fft(signal,[],1);
    recon = complex(zeros(nfrq,NF,NB));
    for n = 1:NB
        for q = 1:nfrq
            recon(q,:,n) = tsig(q,:,n)*A(:,:,q);
        end
    end
    recon = ifft(recon,[],1);
end