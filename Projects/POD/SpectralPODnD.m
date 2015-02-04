function [phi, lambda, f, alpha] = SpectralPODnD(p,FS,cpus)

%Optimized for n spatial dimensions - parallelization is explicit
%Note that for proper behavior, automatic parallel pool creation needs to
%be turned off in the Parallel Computing Toolbox preferences.
%
%Leave FS empty for normalized frequencies
%Leave cpus empty for implicit parallelization
%
%DIM1 is temporal dimension
%DIM2 is spatial dimension
%DIM3 is block dimension (for averaging)

    if ~exist('FS','var')||isempty(FS), FS = 1; end %for normalized frequencies
    
    %Open Matlab Pool if Necessary
    if ~exist('cpus','var')||isempty(cpus), cpus = 0; end
    poolobj = gcp('nocreate');
    if ~isempty(poolobj) && cpus ~= poolobj.NumWorkers
        delete(poolobj);
        poolobj = parpool(cpus);
    elseif isempty(poolobj) && cpus > 0
        poolobj = parpool(cpus);
    end

    %Get frequency information 
    [M,N,NB] = size(p);
    nfreq = ceil((M+1)/2); %number of distinct frequencies
    f = (-nfreq+2:nfreq-1)'*(FS/M);
    
    %Transform into Fourier domain & average
    ph = fft(p,[],1); 

    %Compute eigenvalues, eigenvectors at each frequency & normalize
    phi = complex(zeros(N,N,M)); % organization will be spatial location x mode number x frequency
    lambda = complex(zeros(N,M)); % organization will be mode number x frequency
    alpha = zeros(M,N,NB); % organization will be frequency x mode number x block number
    
    parfor n = 1:M
        %iterate through all blocks, averaging the cross-spectral density
        %(as opposed to the fourier coefficients themselves - you dipshit!)
        R = complex(zeros(N));
        for q = 1:NB
            R = R + ph(n,:,q)'*ph(n,:,q);
        end
        R = (R+R')/(2*NB*M*FS);%Average to remove numerical errors, per Tinney (the matrix is conjugate symmetric), normalize by dt and sample number
        
        %Compute Eigenmodes and eigenvalues
        [phi(:,:,n),eigenvalues] = svd(R);
        lambda(:,n) = diag(eigenvalues);        
    end
    
    for q = 1:NB
        for n = 1:N %we are looping through the POD modes, not the spatial locations
            %Compute temporal POD coefficients
            for k = 1:M
                alpha(k,n,q) =  ph(k,:,q)*squeeze(phi(:,n,k));
            end
        end
    end
    
    %Closes Parallel Pool, if Necessary
    delete(poolobj);
end