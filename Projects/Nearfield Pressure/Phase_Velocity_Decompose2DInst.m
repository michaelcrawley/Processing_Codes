function [subsonic supersonic PSD x f intwaveform data] = Phase_Velocity_Decompose2DInst(src_dir,flist,out_dir,a,windowfun,beta_const,flag)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third (y,x,t). This version only performs the Fourier
    %decomposition along two dimensions, x & t. 
    
    %Constants
    Nx = 22; %number of axial points for upsampling
%     windowfun = @(x) tukeywin(x,0.1);
%     beta_const = 3/sqrt(-log(1e-6));
    D = 0.0254;
    FS = 2e5;
    if ~exist('flag','var'), flag = 'pblocks.p'; end

    for n = 1:length(flist)
        %Load Data
        data = load([src_dir filesep flist{n}]);
        [~,fname] = fileparts(flist{n});
        p = eval(['data.nf.',flag]);
        
        %Get Constants
        [Nt,~,NB] = size(p);
        x = linspace(data.phys.x(1),data.phys.x(end),Nx)*D;
        dx = mean(diff(x));
        dt = 1/FS;
        ds = [dt dx];  
        
        %Initialize Variables
        subsonic = zeros(Nt,Nx,NB);
        supersonic = subsonic;
        intwaveform = zeros(Nt,Nx,NB);
        PSD = subsonic;
        
        for q = 1:NB
            %Interpolate            
            for k = 1:Nt
                intwaveform(k,:,q) = interp1(data.phys.x*D,p(k,:,q),x,'cubic');
            end
        
            %Compute 2-D PSD
            [PSD(:,:,q) f S xm] = PSDN(intwaveform(:,end:-1:1,q),[1 2],ds,windowfun,true); %order of columns needs to be inverted as FFT should be applied as e^(-iwt+ikx)

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
            subsonic(:,end:-1:1,q) = real(iPSDN(S.*subW,[2 1],xm));%order of columns needs to be inverted as iFFT should be applied as e^(iwt-ikx)
            supersonic(:,end:-1:1,q) = real(iPSDN(S.*supW,[2 1],0));
        end
        
        %Save Outputs if necessary 
        if nargout == 0
            save([out_dir filesep fname ' 2Dinst.mat'],'a','windowfun','subsonic','supersonic','data','x','PSD','intwaveform','f');
        end
    end
end