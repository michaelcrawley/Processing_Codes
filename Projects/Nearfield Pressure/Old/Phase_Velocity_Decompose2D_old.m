function [subsonic supersonic angf PSD S vel subW supW] = Phase_Velocity_Decompose2D(waveform,ds,a,windowfun)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third (y,x,t). This version only performs the Fourier
    %decomposition along two dimensions, x & t. 
    
    %Constants
    alpha = 1.0;

    %Compute 3-D PSD
    [PSD f S xm] = PSDN(waveform,[2 3],ds,windowfun,true);
    [M,N,P] = size(S);

    %Change to angular wavenumbers and frequency
    angf = cellfun(@(x) 2*pi*x,f,'uniformoutput',false);

    %Calculate wavenumber magnitude and phase velocity for plane waves
    [kx omega] = meshgrid(angf{1},angf{2});
    vel = abs(omega./kx); %phase velocity
    vel(isnan(vel)) = 0; %get rid of NaN value at kx = 0, w = 0 (will put energy in hydrodynamic portion)
    vel(isinf(vel)) = a; %recapture acoustic energy for kx ~= 0
    vel = repmat(permute(vel,[3 2 1]),[M,1,1]);

    %Determine Weight Vectors (in k,w domain) into subsonically and supersonically traveling waves
    %Decay constant
    dk = mean(diff(f{2}));
    beta = 5*dk;
    epower = 4;
    
    %Create decay about sonic plane
    W = exp(-((vel-a).^epower)/(beta^epower));
       
    %Supersonic weight
    supW = W;
    supW(vel >= alpha*a) = 1;
    
    %Subsonic Weight
    subW = ones(size(supW))-supW;

    %Reconstruct signals
    subsonic = real(iPSDN(S.*subW,[2 3],xm));
    supersonic = real(iPSDN(S.*supW,[2 3],0));
end