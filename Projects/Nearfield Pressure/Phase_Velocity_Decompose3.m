function [subsonic supersonic angf PSD S k vel subW supW] = Phase_Velocity_Decompose3(waveform,ds,a,transform,windowfun)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third.
    %This version is purely for spherical waves
    
    %Compute 3-D PSD
    [PSD f S xm] = PSDN(waveform,[1 2 3],ds,windowfun);
    [M,N,P] = size(S);

    %Change to angular wavenumbers and frequency
    angf = cellfun(@(x) 2*pi*x,f,'uniformoutput',false);

    %Calculate wavenumber magnitude and phase velocity for plane waves
    ky = repmat(angf{1},[1,N]);
    km = repmat(angf{2}.',[M,1]);
    k = transform*[km(:).'; ky(:).']; %transform to orthogonal spatial coordinates
    omega = repmat(permute(angf{3},[2 3 1]),[M N 1]);
    
    %Calculate individual phase velocities (x,y) for spherical waves
    kY = repmat(reshape(k(2,:).',M,N),[1 1 P]);
    kX = repmat(reshape(k(1,:).',M,N),[1 1 P]);
    xvel = omega./kX;
    yvel = omega./kY;
    xvel(isinf(xvel) | isnan(xvel)) = 0;
    yvel(isinf(yvel) | isnan(yvel)) = 0;
    vel = sqrt(xvel.^2+yvel.^2);

    %Determine Weight Vectors (in k,w domain) into subsonically and supersonically traveling waves
    %Decay constant
    df1 = diff(angf{1});
    df2 = diff(angf{2});
    dkv = transform*[df2(1); df1(1)];
    dk = sqrt(dkv.'*dkv);
    beta = 5*dk;
    epower = 4;
    
    %Create decay about sonic plane
    W = exp(-(abs(vel-a).^epower)/(beta^epower));
       
    %Supersonic weight
    supW = W;
    supW(vel >= a) = 1;
    
    %Subsonic Weight
    subW = ones(size(supW))-supW;

    %Reconstruct signals
    subsonic = real(iPSDN(S.*subW,[1 2 3],xm));
    supersonic = real(iPSDN(S.*supW,[1 2 3],0));
end