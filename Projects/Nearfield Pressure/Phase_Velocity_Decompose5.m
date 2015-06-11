function [subsonic supersonic angf PSD S kp ks vel subW supW] = Phase_Velocity_Decompose5(waveform,ds,a,windowfun)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third (y,x,t).
    
    %Constants
    alpha = 1.0;

    %Compute 3-D PSD
    [PSD f S xm] = PSDN(waveform,[1 2 3],ds,windowfun,true);
    [M,N,P] = size(S);

    %Change to angular wavenumbers and frequency
    angf = cellfun(@(x) 2*pi*x,f,'uniformoutput',false);
    
    %Create dimension matrices
    ky = repmat(permute(angf{1}(:),[1 2 3]),[1 N P]);
    kx = repmat(permute(angf{2}(:),[2 1 3]),[M 1 P]);
    omega = repmat(permute(angf{3}(:),[2 3 1]),[M N 1]);    

    %Calculate phase velocity for each dimension
    c_y = omega./ky;
    c_x = omega./kx;
    
    %Calculate total velocity
    vel = sqrt(c_y.^2+c_x.^2);
    vel(isinf(vel)) = a; %recapture acoustic energy for kb ~= 0
    vel(isnan(vel)) = 0; %get rid of NaN value at kb = 0, w = 0 (will put energy in hydrodynamic portion)
    
    %Determine Weight Vectors (in k,w domain) into subsonically and supersonically traveling waves
    %Decay constant
    dk = angf{2}(2)-angf{2}(1); %define based on k_x
    beta = 5*dk;
    epower = 4;
    
    %Create decay about sonic plane
    W = exp(-(abs(vel-a).^epower)/(beta^epower));
       
    %Supersonic weight
    supW = W;
    supW(vel >= alpha*a) = 1;
    
    %Subsonic Weight
    subW = ones(size(supW))-supW;

    %Reconstruct signals
    subsonic = real(iPSDN(S.*subW,[1 2 3],xm));
    supersonic = real(iPSDN(S.*supW,[1 2 3],0));
end