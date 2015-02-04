function Phase_Velocity_Decompose2DInst_v3(src_dir,flist,out_dir,optset)
    %Waveform must have the two spatial dimensions first, and the temporal
    %dimension third (y,x,t). This version only performs the Fourier
    %decomposition along two dimensions, x & t. 
    
    %Set options
    if ~exist('optset','var'), optset = []; end
    if ~isfield(optset,'window'), optset.windowfun = @rectwin; end
    if ~isfield(optset,'flag'), optset.flag = 'pblocks.smp'; end
    if ~isfield(optset,'FS'), optset.FS = 2e5; end
    if ~isfield(optset,'Nx'), optset.Nx = 22; end
    if ~isfield(optset,'PAD'), optset.PAD = [0 4 0]; end
    if ~isfield(optset,'D'), optset.D = 0.0254; end
    if ~isfield(optset,'temp_flag'), optset.temp_flag = 'amb'; end
    if ~isfield(optset,'beta_const'), optset.beta_const = 1e-5; end
    dt = 1/optset.FS;
    matversion = 1.2;
    fields = [1 2 3 7 14];

    for k = 1:length(flist)
        %load file
        tmpdata = load([src_dir filesep flist{k}]);
        [~,fname] = fileparts(flist{k});
        sig = eval(['tmpdata.nf.',optset.flag]);
        if strcmpi(optset.temp_flag,'amb')
            a = tmpdata.phys.a;
        elseif strcmpi(optset.temp_flag,'jet')
            Temp = tmpdata.phys.To/(1+0.2*tmpdata.phys.M^2);
            a = sqrt(1.4*287*Temp);
        end
               
        %Get Parameters
        N = size(sig);
        x = tmpdata.phys.x*optset.D; %convert to meters
        int_x = linspace(x(1),x(end),optset.Nx); %interpolate onto regular grid
        dx = mean(diff(int_x));              
        L = 2.^(nextpow2(N)+optset.PAD); %Find zero-pad parameters
                
        %Grab any additional helpful data
        ff = eval(['tmpdata.ff.',optset.flag]);
        
        %Initialize large matrices
        intwaveform = zeros(N(1),optset.Nx,N(3));
        subsonic = intwaveform;
        supersonic = intwaveform;
        
        %Create window (time only)
        wndo_t = repmat(optset.windowfun(N(1)),[1 optset.Nx]);
        
        for q = 1:N(3)
            %Interpolate            
            tmp = zeros(N(1),optset.Nx);
            for n = 1:N(1)
                tmp(n,:) = interp1(x,sig(n,:,q),int_x,'pchip');
            end
            intwaveform(:,:,q) = tmp;
            
            %Apply window (in time only) , zero-pad
            tmp = intwaveform(:,:,q).*wndo_t;
            swap = zeros(L(1),L(2));
            swap(1:N(1),1:optset.Nx) = tmp;
            tmp = swap; 
        
            %Compute 2-D PSD
            [PSD(:,:,q), f, S, xm] = PSDN(tmp(:,end:-1:1),[1 2],[dt dx],@rectwin,true); %order of columns needs to be inverted as FFT should be applied as e^(-iwt+ikx)

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
                beta = optset.beta_const*dk;

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
            sub = sub(1:N(1),1:optset.Nx);
            sup = sup(1:N(1),1:optset.Nx); 
            
            %Get rid of time-windowing effects;
            subsonic(:,:,q) = sub;
            supersonic(:,:,q) = sup;
        end
        
        
        PSD = mean(PSD,3);
        
        %Organize variables
        filter_params = optset;
        filter_params.a = a;
        filter_params.x = int_x;
        filter_params.file = mfilename;
        
        data.ff = ff;
        data.intwaveform = intwaveform;
        data.subsonic = subsonic;
        data.supersonic = supersonic;
        data.PSD = PSD;
        
        daq_params.phys = tmpdata.phys;
        daq_params.phys.D = optset.D;
        daq_params.filename = tmpdata.filename;
        daq_params.trigger = tmpdata.trigger;
        daq_params.phys.FS = optset.FS;
        
        %Get new filename
        [~,fname] = fileparts(flist{k});
        parts = regexp(fname,'_','split'); 
        short = parts(fields);
        newname =[sprintf('%s_',short{1:end-1}),short{end},' 2DFourierInst.mat'];
        
        %Save Output
        save([out_dir,filesep,newname],'data','filter_params','daq_params','matversion');
    end
end