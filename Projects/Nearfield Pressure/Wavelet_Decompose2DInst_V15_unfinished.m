function Wavelet_Decompose2DInst_V15(src_dir,flist,out_dir,a,windowfun,flag,cpus)

    %Constants
    if ~exist('windowfun','var')||isempty(windowfun), windowfun = @tukeywin; end
    if ~exist('flag','var')||isempty(flag), flag = 'pblocks.p'; end
    if ~exist('cpus','var')||isempty(cpus), cpus = 1; end
    
    %Parameters
    mother = 'morlet';
    param = [6 6];
    dj = [1/3 1/3];
    dt = 1/2e5;
    D = 0.0254;
    Nx = 22;
    DownFS =  8; %Frequency downsampling -> FS/DownFS
    dt = dt*DownFS;
    DownT = 8; %Period downsampling N = 8192 -> N = 4096
    PAD = [0 0 0];
    
    %Open Matlab Pool if Necessary
    if cpus > 1
        if matlabpool('size') > 0
            matlabpool close;
        end
        matlabpool(cpus);
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
        L = 2.^(nextpow2(N)+PAD); %Find zero-pad parameters
        
        
        %Grab any additional helpful data
        ff = eval(['data.ff.',flag]);
        ff = ff(1:DownFS:end/DownT,:,:);
        ff = ff(t_indx,:,:); %cut down far-field spectra to match output of decomposition
        
        %Initialize large matrices
        intwaveform = zeros(length(t_indx),Nx,N(3));
        subsonic = intwaveform;
        supersonic = intwaveform;
        
        if k == 1
            %Compute Analyzing Wavelets - scales{1} are temporal, scales{2}
            %are spatial
            [sub_filt,sup_filt,C_delta,total_filt] = DaughterWavelets2(L,[dt dx],dj,a,mother,param);
        end
        
        %Compute Decompositions
        disp(['Processing File: ',fname,'...']);
        tic;
        for n = 1:N(3)            
            %Interpolate
            tmp = interp2(X,T,sig(:,:,n),new_X,new_T);
            intwaveform(:,:,n) = tmp(t_indx,:);
            
            %Apply window, zero-pad
            wndo_t = repmat(windowfun(N(1)),[1 Nx]); %window along t-dimension (DIM1)
            tmp = tmp.*wndo_t;
            swap = zeros(L(1),L(2));
            swap(1:N(1),1:Nx) = tmp;
            tmp = swap; clear swap;         
            
            %Filtered Reconstructions
            xh = fftn(tmp);
            tmp = ifftn(xh.*conj(sub_filt))/C_delta;
            subsonic(:,:,n) = real(tmp(t_indx,1:Nx));
            tmp = ifftn(xh.*conj(sup_filt))/C_delta;
            supersonic(:,:,n) = real(tmp(t_indx,1:Nx));
            tmp = ifftn(xh.*conj(total_filt))/C_delta;
            full(:,:,n) = real(tmp(t_indx,1:Nx));
        end
        %Update User
        disp(['Complete...',' Processing time: ',num2str(round(toc)),' seconds...']);
        
        %Save Output
        save([out_dir,filesep,fname,' 2DWave TEST3.mat'],'intwaveform','subsonic','supersonic','a','windowfun','int_x','data','ff','full','sub_filt','sup_filt','total_filt');
    end
    
    %Closes Parallel Pool, if Necessary
    if cpus > 1
        matlabpool('close');
    end
end

function [sub_filt,sup_filt,C_delta,total] = DaughterWavelets2(L,dt,dj,ca,mother,param)
    %Grab Mother Wavelet info
    if ~exist('dj','var')||isempty(dj), dj = [0.125 0.125]; end
    if ~exist('mother','var')||isempty(mother), mother = 'morlet'; end
    if ~exist('param','var')||isempty(param), param = [6 6]; end
    if ~exist('ca','var')||isempty(ca), flag = false; else flag = true; end
    [~,fcoef] = MotherWavelets('ST2',mother,param);
    
    %Compute scales
    c0 = 2*dt(2)*fcoef(2)/(L(1)*dt(1)*fcoef(1)); %initial velocity scales
    cN = ceil((log2(L(1)*L(2))-2)/dj(1)); %number of velocity scales
    scales{1} = c0*2.^(dj(1)*(0:cN)); %velocity scales
    s0 = (1/2/pi)*sqrt(prod(fcoef))*sqrt(prod(dt)); %initial spatial scales
    sN = ceil((log2(sqrt(L(1)*L(2)))-1)/dj(2)); %number of spatial scales
    scales{2} = s0*2.^(dj(2)*(0:sN)); %spatial scales
    
    %Recompute velocity scales to hit acoustic velocity exactly
    if flag
        X = round(cN*log2(ca/scales{1}(1))/log2(scales{1}(end)/scales{1}(1)));
        dj2 = log2(ca/scales{1}(1))/X;
        scales{1} = c0*2.^(dj2*(0:cN));
    end
    sub_indx = find(scales{1}+1e3*eps < ca,1,'last'); %get scale(velocity) index for subsonic velocities
    
    %compute normalized FFT axes
    omega_t = 2*pi*(1:floor(L(1)/2))/L(1)/dt(1);
    omega_x = 2*pi*(1:floor(L(2)/2))/L(2)/dt(2);
    omegaks{1} = [0,omega_t,-fliplr(omega_t(1:floor((L(1)-1)/2)))];
    omegaks{2} = [0,omega_x,-fliplr(omega_x(1:floor((L(2)-1)/2)))];
    [K,T] = meshgrid(omegaks{2},omegaks{1}); %order is flipped because meshgrid assumes input is DIM2,DIM1
    
    %Compute all daughter wavelets
    cN1 = cN+1;
    sN1 = sN+1;
    
    wavelet = MotherWavelets('ST2',mother,param);
    %Calculate subsonic-filtered map - method 1
    sub_filt1 = zeros(L(1),L(2));
    tic;
    disp('Computing Subsonic Map....');
    for j = 1:sub_indx
        for k = 1:sN1
            sub_filt1 = sub_filt1 + wavelet(scales{2}(k)*T/sqrt(scales{1}(j)),scales{2}(k)*K*sqrt(scales{1}(j)))/sqrt(scales{1}(j)*scales{2}(k));
        end
    end
    disp(['Complete...',' Processing time: ',num2str(round(toc)),' seconds...']);
    
    %calculate subsonic-filtered map - method 2
    sub_filt2 = zeros(L(1),L(2),length(scales{1}(1:sub_indx)),length(scales{2}));
    tic;
    disp('Computing Subsonic Map....');
    for j = 1:sub_indx
        for k = 1:sN1
            sub_filt2(:,:,j,k) = wavelet(scales{2}(k)*T/sqrt(scales{1}(j)),scales{2}(k)*K*sqrt(scales{1}(j)))/(scales{1}(j)*scales{2}(k));
        end
    end
    sub_filt2 = trapz(scales{2},sub_filt2,4);
    sub_filt2 = trapz(scales{1}(1:sub_indx),sub_filt2,3);
    disp(['Complete...',' Processing time: ',num2str(round(toc)),' seconds...']);
    
    %Calculate sonic/supersonic-filtered map
    sup_filt = zeros(L(1),L(2));
    tic;
    disp('Computing Supersonic Map...');
    for j = sub_indx+1:cN1
        for k = 1:sN1
            wavelet = MotherWavelets('ST2',mother,param);
            sup_filt = sup_filt + wavelet(scales{2}(k)*T/sqrt(scales{1}(j)),scales{2}(k)*K*sqrt(scales{1}(j)))/sqrt(scales{1}(j)*scales{2}(k));
        end
    end
    disp(['Finished...Compute Time: ',num2str(toc)]);
    
    %Calculate total energy
    total = sub_filt + sup_filt;
    C_delta = mean(real(total(:)));
end


function [wavelet,fcoef] = MotherWavelets(dim,mother,m)
%This function returns a function handle, 'wavelet', given a dimension and
%string for the name of the mother wavelet. A parameter to modify the base
%mother wavelet can also be provided. For standard mother wavelets, the
%dimension provided will be an integer, for Spatio-Temporal mother
%wavelets, the dimension will be a string beginning with 'ST' and
%including the total number of dimensions.

    if ~exist('m','var'), m = []; end

    switch dim
        case 1
            switch lower(mother)
                case 'morlet'
                    if isempty(m), m = 6; end
                    fcoef = (4*pi)/(m + sqrt(2 + m^2));
                    wavelet = @(k) (pi^-0.25).*exp(-0.5*((k - m).^2)).*(k > 0); 
                case 'paul'
                    if isempty(m), m = 4; end
                    fcoef = 4*pi/(2*m+1);
                    wavelet = @(k) (2^m/sqrt(m*factorial(2*m-1))).*(k.^m).*exp(-(k).*(k > 0)).*(k > 0);
                case 'dog'
                    if isempty(m), m = 2; end
                    fcoef = 2*pi/sqrt(m+0.5);
                    wavelet = @(k) -(1i^m)/sqrt(gamma(m+0.5))*(k.^m).*exp(-0.5*k.^2);
                otherwise
                    error('Undefined Mother Wavelet');
            end

        case 'ST2'
            switch lower(mother)
                case 'morlet'
                    if isempty(m), m = [6 6]; end
                    fcoef = m + sqrt(2 + m.^2);
                    wavelet = @(omega,k) (pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k > 0) + ... %positive wavenumbers
                                            (pi^-0.5)*exp(-0.5*((k + m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k < 0); %negative wavenumbers 
                otherwise
                    error('Undefined Mother Wavelet');
            end
        otherwise
            error('Incorrect Dimension Definition');
    end
end