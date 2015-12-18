function [recon,A,f] = MultiSpectralLSE(field,signal,I,FS)

%This is a highly-specialized code, designed specifically to upsample
%non-time resolved PIV data (conditional field) using a generic set of
%microphone data (unconditional signal). This code utilizes the algorithm
%laid out by Tinney 2008 JFM.
%Inputs:
%           field:  Conditional field to be reconstructed, 3-D matrix where 
%                   DIM1 corresponds to sparse time (i.e. only a single 
%                   time index is populated), DIM2 corresponds to spatial
%                   locations, and DIM3 corresponds to blocks.
%           signal: Unconditional reference signal, 3-D full matrix.
%           FS:     Sampling frequency (optional)
%Outputs:
%           recon:  Conditional reconstruction of field
%           A:      Frequency-dependent reconstruction matrix
%           f:      Frequency axis 

    if ~exist('FS','var')||isempty(FS), FS = 1; end
    
    [BS,NS,NB] = size(signal);
    NF = size(field,2);
    
    %Verify that field is temporally sparse
%     [I,~,~] = find(field); %must use all outputs, so that I only corresponds to first dim
%     I = unique(I);
%     if length(I) > 1
%         error('Field is not temporally-sparse!');
%     end
    
    unconditional = zeros(NS,NS,BS);
    conditional = zeros(NS,NF,BS);    
    for n = 1:NB
        for q = 1:NS
            for k = 1:NS
                unconditional(q,k,:) = unconditional(q,k,:) + permute(fft(flipud(signal(:,q,n))*signal(I,k,n),[],3)/NB,[2 3 1]); %remember that the signal gets flipped in the convolution
            end
            for k = 1:NF
                conditional(q,k,:) = conditional(q,k,:) + permute(fft(flipud(signal(:,q,n))*field(n,k),[],3)/NB,[2 3 1]);
            end
        end        
    end

    %Transform into Fourier domain and compute estimation coefficients at
    %each frequency
%     unconditional = fft(unconditional,[],3);
%     unconditional = mean(unconditional,4);
%     conditional = fft(conditional,[],3);
%     conditional = mean(conditional,4);
    if mod(BS,2) %signal length is odd
        Nc = (BS-1)/2;
        f = ifftshift(-Nc:Nc)*FS/BS;
    else
        Nc = BS/2;
        f = [0,1:Nc,-(Nc-1):1:-1]*FS/BS;
    end
    
    A = complex(zeros(NS,NF,BS));
    for n = 1:BS
        A(:,:,n) = unconditional(:,:,n)\conditional(:,:,n);
    end
    A(isnan(A)|isinf(A)) = 0;%fix for badly conditioned sets
    clear conditional unconditional;
    
    %Compute reconstruction in Fourier domain, transform into time domain
    tsig = fft(signal,[],1);
%     recon = complex(zeros(BS,NF,NB));
    recon = complex(zeros(BS,NF));
    for n = 1:1 %NB
        for q = 1:BS
            recon(q,:,n) = tsig(q,:,n)*A(:,:,q);
        end
    end
    recon = ifft(recon,[],1);
end