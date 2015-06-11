function [subsonic supersonic angf PSD S k vel subW supW] = Phase_Velocity_Decompose2(waveform,ds,a,transform,windowfun)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third. Unlike the first iteration, this code uses the cutoff
    %defined by Cabana et al 2008: kx^2+ky^2 < (omega/a)^2

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
    kb = sqrt(k(1,:).^2+k(2,:).^2); %get magnitude of wavenumber vectors
    kb = reshape(kb,M,N);
    kb = repmat(kb,[1 1 P]);
    vel = omega./kb; %phase velocity
    vel(isinf(vel) | isnan(vel)) = 0; %get rid of inf value at kb = 0
    kY = repmat(reshape(k(2,:).',M,N),[1 1 P]);
    kX = repmat(reshape(k(1,:).',M,N),[1 1 P]);
    ineq = (omega/a).^2;  %kx^2 + ky^2 < ineq for sonic/supersonic
    
    %Determine Weight Vectors (in k,w domain) into subsonically and supersonically traveling waves
    %Decay constant
    df1 = diff(f{1});
    df2 = diff(f{2});
    dkv = transform*[df2(1); df1(1)];
    dk = sqrt(dkv.'*dkv);
    beta = 5*dk;
    epower = 4;
    
    %Create decay about sonic plane
    W = exp(-(sqrt(kX.^2+kY.^2-ineq).^epower)/(beta^epower));
       
    %Supersonic weight
    supW = W;
    flag = kX.^2+kY.^2 <= ineq;
    supW(flag) = 1;
    
    %Subsonic Weight
    subW = ones(size(supW))-supW;

    %Reconstruct signals
    subsonic = real(iPSDN(S.*subW,[1 2 3],xm));
    supersonic = real(iPSDN(S.*supW,[1 2 3],0));
end