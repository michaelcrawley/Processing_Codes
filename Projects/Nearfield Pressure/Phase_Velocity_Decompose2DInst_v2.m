function [subsonic supersonic PSD x f intwaveform data] = Phase_Velocity_Decompose2DInst_v2(src_dir,flist,out_dir,a,flag)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third (y,x,t). This version only performs the Fourier
    %decomposition along two dimensions, x & t. 
    
    %Constants
    Nx = 22; %number of axial points for upsampling
    windowfun = @rectwin;
    beta_const = (1e-5);
    D = 0.0254;
    DownFS = 1;
    DownT = 1;
    FS = 2e5/DownFS;
    XPAD = 5;
    if ~exist('flag','var'), flag = 'pblocks.smp'; end
    matversion = 1.1;

    for n = 1:length(flist)
        %Load Data
        data = load([src_dir filesep flist{n}]);
        [~,fname] = fileparts(flist{n});
        p = eval(['data.nf.',flag]);
        
        %Downsample
        p = p(1:DownFS:end/DownT,:,:);
        
        %Get Constants
        [Nt,~,NB] = size(p);
        x = linspace(data.phys.x(1),data.phys.x(end),Nx)*D;
        dx = mean(diff(x));
        dt = 1/FS;
        ds = [dt dx]; 
        t = (0:Nt-1)*dt;
        x = data.phys.x*D; %convert to meters
        [X,T] = meshgrid(x,t); %for interpolating
        int_x = linspace(x(1),x(end),Nx); %interpolate onto regular grid
        [new_X,new_T] = meshgrid(int_x,t);  
        
        %Create window (time only)
        wndo_t = repmat(windowfun(Nt),[1 Nx]);
        t_indx = find(abs(windowfun(Nt) - 1) < 10*eps);
        
        %Initialize Variables
        subsonic = zeros(length(t_indx),Nx,NB);
        supersonic = subsonic;
        intwaveform = subsonic;
        
        for q = 1:NB
            %Interpolate            
            for k = 1:Nt
                intwaveform(k,:,q) = interp1(data.phys.x*D,p(k,:,q),int_x,'pchip');
            end
            
            %Apply window (in time only) 
            sig = intwaveform(:,:,q).*wndo_t;            
            
            %Zero-Pad
            L1 = 2^(nextpow2(Nt)); %pad in time
            L2 = 2^(nextpow2(Nx)+XPAD); %increase resolution in x-direction only
            swap = zeros(L1,L2);
            swap(1:Nt,1:Nx) = sig;
            sig = swap; clear swap;
        
            %Compute 2-D PSD
            [PSD(:,:,q) f S xm] = PSDN(sig(:,end:-1:1),[1 2],ds,@rectwin,true); %order of columns needs to be inverted as FFT should be applied as e^(-iwt+ikx)

            if q == 1 %These values will be identical for all blocks
                %Change to angular wavenumbers and frequency
                angf = cellfun(@(x) 2*pi*x,f,'uniformoutput',false);

                %Calculate wavenumber magnitude and phase velocity for plane waves
                [kx omega] = meshgrid(angf{2},angf{1});
                vel = abs(omega./kx); %phase velocity
                vel(isnan(vel)) = 0; %get rid of NaN value at kx = 0, w = 0 (will put energy in hydrodynamic portion)
                vel(isinf(vel)) = a; %recapture acoustic energy for kx ~= 0
                ka = omega/a; %acoustic wavenumber
                dka = (abs(kx)-abs(ka))/2/pi;

                %Determine Weight Vectors (in k,w domain) into subsonically and supersonically traveling waves
                %Decay constant
                dk = f{2}(2);
                beta = beta_const*dk;

                %Create decay about sonic plane
                W = exp(-(dka/beta).^2);

                %Supersonic weight
                supW = W;
                supW(vel >= a) = 1;

                %Subsonic Weight
                subW = W;
                subW(vel <= a) = 1;
%                 subW = ones(size(supW))-supW;
            end

            %Reconstruct signals
            sub = fliplr(real(iPSDN(S.*subW,[2 1],xm)));%order of columns needs to be inverted as iFFT should be applied as e^(iwt-ikx)
            sup = fliplr(real(iPSDN(S.*supW,[2 1],0)));
            
            %Get rid of zero-padded cells
            sub = sub(1:Nt,1:Nx);
            sup = sup(1:Nt,1:Nx);  
            
            %Get rid of time-windowing effects;
            subsonic(:,:,q) = sub(t_indx,:);
            supersonic(:,:,q) = sup(t_indx,:);
        end
        
        PSD = mean(PSD,3);
        intwaveform = intwaveform(t_indx,:,:);
        ff = eval(['data.ff.',flag,'(t_indx,:,:)']);
        
        %Save Outputs if necessary 
        if nargout == 0
            save([out_dir filesep fname ' 2Dinst.mat'],'a','windowfun','subsonic','supersonic','data','x','PSD','intwaveform','f','ff','flag');
        end
    end
end