function [phi,psi] = DirectSpectralEstimation(p,x,y,a,FS)
%Based off of the method proposed by Papamoschou in JSV 2011
%Inputs:
%           p:          signal (blocksize x channels x number of blocks)
%           x:          axial location of microphones (m)
%           y:          radial location of microphones (m)
%           a:          ambient speed of sound (m/s)
%           FS:         sampling frequency (Hz)
%Outputs:
%           phi:        coherence-based noise-source distribution
%           psi:        spectrum-based noise-source distribution

    %constants
    qmax = 100; %max number of iterations for minimization

    %Get processing info
    N = size(p);
    J = N(2)^2 - N(2) + 1; %number of indepedent values in complex-coherence matrix
    K = J; %number of assumed noise sources
    nfrq = ceil((N(1)+1)/2);
    f = 2*pi*(0:nfrq-1)*(FS/N(1)); %angular frequency
    
    %Define noise source locations
    xs = linspace(0,x(end),K); %defined as axial array, on jet axis
    [Xs,X] = meshgrid(xs,x);
    Y = repmat(y(:),[1,K]);
    l = sqrt((X-Xs).^2 + Y.^2); %therefore, the distance from noise source k to microphone m is located at l(m,k)
    tau = l/a;
    
    %initialize variables
    phi = zeros(nfrq-1,K);
    psi = zeros(nfrq-1,K,N(2));
    
    %Cycle through frequencies
    ph = fft(p,[],1);    
    for n = 2:nfrq %skip zero-frequency (first index of fft)        
        %Compute complex-coherence
        tmp = squeeze(ph(n,:,:));
        cpsd = tmp*tmp'/N(3); %note that we want to use conjugate-transpose here!
        psd = diag(cpsd);
        tmp = repmat(sqrt(psd),[1,N(2)]);
        norm = tmp*tmp'; %tmp is real (by definition), so the conjugate-transpose operator doesn't matter
        coh = cpsd./norm; %complex-coherence, so negatives are possible!
        
        %Compute array response matrix
        Tm = exp^(1i*f(n)*tau);
        
        %Iterate to find 'true' distribution for a                
        a = ones(length(xs),1);%initial condition for noise source distribution
        for q = 1:qmax
            
        end
        
        %Compute source distributions
        phi(n-1,:) = a.^2; %coherence based noise-source distribution (NEED TO GET RID OF DISCRETIZATION!!! See Eq. 28 - 29)
    end

end