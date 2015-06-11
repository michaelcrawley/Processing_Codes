function [subsonic supersonic angf PSD S kp vel subW supW] = Phase_Velocity_Decompose4(waveform,ds,a,transformp,transforms,windowfun)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third (y,x,t).
    
    %Constants
    alpha = 1.0;

    %Compute 3-D PSD
    [PSD f S xm] = PSDN(waveform,[1 2 3],ds,windowfun,true);
    [M,N,P] = size(S);

    %Change to angular wavenumbers and frequency
    angf = cellfun(@(x) 2*pi*x,f,'uniformoutput',false);

    %Calculate wavenumber magnitude and phase velocity for plane waves
    [km ky] = meshgrid(angf{2},angf{1});
%     ky = repmat(angf{1},[1,N]);
%     km = repmat(angf{2}.',[M,1]);
    kp = transformp*[km(:).'; ky(:).']; %transform to orthogonal spatial coordinates for plane waves
    kb = sqrt(kp(1,:).^2+kp(2,:).^2); %get magnitude of wavenumber vectors
    kb = reshape(kb,M,N);
    kb = repmat(kb,[1 1 P]);
    omega = repmat(permute(angf{3},[2 3 1]),[M N 1]);
    vel = abs(omega./kb); %phase velocity
    vel(isnan(vel)) = 0; %get rid of NaN value at kb = 0, w = 0 (will put energy in hydrodynamic portion)
    
    %Calculate individual phase velocities (x,y) for spherical waves
    ks = transforms*[km(:).'; ky(:).']; %transform to orthogonal spatial coordinates for spherical waves
    kY = repmat(reshape(ks(2,:).',M,N),[1 1 P]);
    kX = repmat(reshape(ks(1,:).',M,N),[1 1 P]);
    xvel = abs(omega./kX);
    yvel = abs(omega./kY);
    vel((xvel >= a) & (yvel >= a)) = a; %if both the x and y phase velocities are sonic/supersonic, the wave must be curved and acoustic
    vel(isinf(vel)) = a; %recapture acoustic energy for kb ~= 0

    %Determine Weight Vectors (in k,w domain) into subsonically and supersonically traveling waves
    %Decay constant
    df1 = diff(f{1});
    df2 = diff(f{2});
    dkv = transforms*[df2(1); df1(1)];
    dk = sqrt(dkv.'*dkv);
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
    subsonic = real(iPSDN(S.*subW,[1 2 3],xm));
    supersonic = real(iPSDN(S.*supW,[1 2 3],0));
end