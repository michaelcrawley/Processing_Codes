function [phi lambda e_lambda e_phi PODrec PSD] = SpectralPOD1D(p)

%For 1 spatial dimension only
%DIM1 is temporal dimension
%DIM2 is spatial dimension
%DIM3 is block dimension (for averaging)

    %Transform into Fourier domain & average
    ph = fft(p,[],1);
    ph = mean(ph,3);
    
    %Get frequency information & cut out negative frequencies
    [M,N] = size(ph);
    nfreq = ceil((M+1)/2);
    ph = ph(1:nfreq,:);
    
    %Compute cross-correlation matrix    
    R = zeros(N,N,nfreq);
    for n = 1:N
        R(n,:,:) = permute(repmat(conj(ph(:,n)),1,N).*ph,[2 1]); %note that we need the complex conjugate here
    end
    
    %Compute eigenvalues, eigenvectors at each frequency & normalize
    for n = 1:nfreq
        R(:,:,n) = (R(:,:,n)+R(:,:,n)')/2;%Average to remove numerical errors, per Tinney (the matrix is conjugate symmetric)
        [phi(:,:,n),eigenvalues] = svd(R(:,:,n));
        lambda(:,n) = diag(eigenvalues)/N;
        phi(:,:,n) = sqrt(repmat(lambda(:,n).',N,1)).*phi(:,:,n);
    end
    
    %Check reconstruction
    for n = 1:nfreq
        tmp = N*abs(phi(:,:,n)).^2;
        PODrec(:,n) = abs(sum(tmp,2));
    end
       
    %Normalize eigenvalues
    lambda = real(lambda); %take real part only, because Ching-Wen says so
    
    %Verify 
    for n = 1:nfreq
        e_lambda(n) = sum(lambda(:,n));
        e_phi(n) = sum(diag(R(:,:,n)))/N;
        PSD(n,:) = diag(R(:,:,n));
    end
end