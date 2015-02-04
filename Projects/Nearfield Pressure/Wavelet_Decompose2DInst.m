function Wavelet_Decompose2DInst(src_dir,flist,out_dir,a,windowfun,flag,cpus)

    %Constants
    if ~exist('windowfun','var')||isempty(windowfun), windowfun = @tukeywin; end
    if ~exist('flag','var')||isempty(flag), flag = 'pblocks.p'; end
    if ~exist('cpus','var')||isempty(cpus), cpus = 1; end
    
    %Parameters
    mother = 'morlet';
    param = [6 6];
    dj = [1/13 1/11];
    dt = 1/2e5;
    D = 0.0254;
    Nx = 22;
    DownFS =  2; %Frequency downsampling -> FS/DownFS
    dt = dt*DownFS;
    DownT = 1; %Period downsampling N = 8192 -> N = 4096
    
    %Open Matlab Pool if Necessary
    if cpus > 1
        cp = gcp('nocreate');
        if ~isempty(cp)
            delete(cp);
        end
        cp = parpool(cpus);
        
    end
    
    for k = 1:length(flist)
        %load file
        data = load([src_dir filesep flist{k}]);
        [~,fname] = fileparts(flist{k});
        sig = eval(['data.nf.',flag]);
        
        %Downsample
        sig = sig(1:DownFS:end/DownT,:,:);
        
        %Get Parameters
        N = size(sig);
        t_indx = find(abs(windowfun(N(1)) - 1) < 10*eps); %cut down to account for gaussian in tukey window
        t = (0:N(1)-1)*dt;
        x = data.phys.x*D; %convert to meters
        [X,T] = meshgrid(x,t); %for interpolating
        int_x = linspace(x(1),x(end),Nx); %interpolate onto regular grid
        dx = mean(diff(int_x));
        [new_X,new_T] = meshgrid(int_x,t);
        
        %Grab any additional helpful data
        ff = eval(['data.ff.',flag]);
        ff = ff(1:DownFS:end/DownT,:,:);
        ff = ff(t_indx,:,:); %cut down far-field spectra to match output of decomposition
        
        %Initialize large matrices
        intwaveform = zeros(length(t_indx),Nx,N(3));
        subsonic = intwaveform;
        supersonic = intwaveform;
        
        %Compute Decompositions
        disp(['Processing File: ',fname,'...']);
        for n = 1:N(3)
            tic;
            %Interpolate
            tmp = interp2(X,T,sig(:,:,n),new_X,new_T);
            intwaveform(:,:,n) = tmp(t_indx,:);
            
            %Forward Wavelet Transform            
            [wave,scales,omegaks] = STCWT1D(tmp,[dt dx],dj,windowfun,a,mother,param,'-real');
            
            %Following computations only need to be done once
            if n == 1
                %Calculate Normalization Constant, C_delta, for inverse
                %transform
                C_delta = STCWT1D_Cdelta(scales,omegaks,mother,param);

                %find Index for last subsonic phase velocities
                I = find(scales{1} < a,1,'last');
                
                %Cut cell arrays for reconstruction                
                sub_scales{1} = scales{1}(1:I); %grab subsonic phase velocities only
                sub_scales{2} = scales{2};  
                
                sup_scales{1} = scales{1}(I+1:end); %grab supersonic phase velocities only
                sup_scales{2} = scales{2};                
            end

            %Reconstructions
            subsonic(:,:,n) = iSTCWT1D(wave(1:I,:,t_indx,:),sub_scales,[],[],[],C_delta);
            supersonic(:,:,n) = iSTCWT1D(wave(I+1:end,:,t_indx,:),sup_scales,[],[],[],C_delta);
            
            %Update User
            disp([num2str(round(100*n/N(3))),'% Complete...',' Processing time: ',num2str(round(toc)),' seconds...']);
        end
        
        %Save Output
        save([out_dir,filesep,fname,' 2DWave.mat'],'intwaveform','subsonic','supersonic','a','windowfun','int_x','data','ff');
    end
    
    %Closes Parallel Pool, if Necessary
    if cpus > 1
        delete(cp);
    end
end