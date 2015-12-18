function Wavelet_Decompose2DInst_V4_2(src_dir,flist,out_dir,a,flag,cpus)
%Computes STCWT in the physical domain, whereas V3 computed in the Fourier
%domain! In this code, the wavelet filter is first computed in the Fourier
%domain, and then transformed into the physical domain.
    tic;
    %Constants
    if ~exist('flag','var')||isempty(flag), flag = 'pblocks.p'; end
    if ~exist('cpus','var')||isempty(cpus), cpus = 1; end
    
    %Parameters
    mother = 'morlet';
    param = [6 6];
    dj = [1/18 1/18];
    dt = 1/2e5;
    D = 0.0254;
    Nx = 21*2+1;
    DownFS =  1; %Frequency downsampling -> FS/DownFS
    dt = dt*DownFS;
    DownT = 1; %Period downsampling N = 8192 -> N = 4096
    PAD = [0 5 0];
    
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
        
        %Initialize large matrices
        intwaveform = zeros(N(1),Nx,N(3));
        subsonic = intwaveform;
        supersonic = intwaveform;
        
        if k == 1
            %Compute Analyzing Wavelets - scales{1} are temporal, scales{2}
            %are spatial
            [sub_filt,sup_filt,C_delta] = DaughterWavelets(L,[dt dx],dj,a,mother,param);
            
            %Cut out padding
            cut = floor((L(1:2)-[N(1),Nx])/2);
            sub_filt = sub_filt(cut(1)+1:(end-cut(1)),cut(2)+1:(end-cut(2)));
            sup_filt = sup_filt(cut(1)+1:(end-cut(1)),cut(2)+1:(end-cut(2)));
        end
        
        %Compute Decompositions
        disp(['Processing File: ',fname,'...']);
        for n = 1:N(3)            
            %Interpolate
            intwaveform(:,:,n) = interp2(X,T,sig(:,:,n),new_X,new_T,'spline');
            
            %Filtered Reconstructions 
            subsonic(:,:,n) = conv2(intwaveform(:,:,n),sub_filt,'same')/C_delta;
            supersonic(:,:,n) = conv2(intwaveform(:,:,n),sup_filt,'same')/C_delta;
        end
        
        %Save Output
        save([out_dir,filesep,fname,' 2DWaveDecomp.mat'],'intwaveform','subsonic','supersonic','a','int_x','data','ff','sub_filt','sup_filt');
    end
    
    %Closes Parallel Pool, if Necessary
    if cpus > 1
        matlabpool('close');
    end
    disp(['Complete...',' Processing time: ',num2str(round(toc)),' seconds...']);
end

function [sub_filt,sup_filt,C_delta,total] = DaughterWavelets(L,dt,dj,ca,mother,param)
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
    [C,S] = ndgrid(scales{1},scales{2}); %order is flipped because meshgrid assumes input is DIM2,DIM1
    
    %Compute all daughter wavelets

    wavelet = MotherWavelets('ST2',mother,param);
    
    %Calculate subsonic-filtered map
    sub_filt = cell(L(1),L(2));
    sup_filt = cell(L(1),L(2));
    disp('Computing Phase Maps....');
%     for j = 1:length(omegaks{1})
%         for k = 1:length(omegaks{2})
%             %Subsonic map
%             tmp = wavelet(S(1:sub_indx,:)*omegaks{1}(j)./sqrt(C(1:sub_indx,:)),S(1:sub_indx,:)*omegaks{2}(k).*sqrt(C(1:sub_indx,:)))./(C(1:sub_indx,:).*S(1:sub_indx,:));
%             sub_filt(j,k) = trapz(scales{1}(1:sub_indx),trapz(scales{2},tmp,2));
%             
%             %Supersonic map
%             tmp = wavelet(S(sub_indx+1:end,:)*omegaks{1}(j)./sqrt(C(sub_indx+1:end,:)),S(sub_indx+1:end,:)*omegaks{2}(k).*sqrt(C(sub_indx+1:end,:)))./(C(sub_indx+1:end,:).*S(sub_indx+1:end,:));
%             sup_filt(j,k) = trapz(scales{1}(sub_indx+1:end),trapz(scales{2},tmp,2));
%         end
%     end

    Nmax = length(omegaks{1})*length(omegaks{2});
    parfor q = 1:Nmax
            
            [j,k] = ind2sub( [length(omegaks{1}) length(omegaks{2})], q );
        
            %Subsonic map
            tmp = wavelet(S(1:sub_indx,:)*omegaks{1}(j)./sqrt(C(1:sub_indx,:)),S(1:sub_indx,:)*omegaks{2}(k).*sqrt(C(1:sub_indx,:)))./(C(1:sub_indx,:).*S(1:sub_indx,:));
            sub_filt{q} = trapz(scales{1}(1:sub_indx),trapz(scales{2},tmp,2));
            
            %Supersonic map
            tmp = wavelet(S(sub_indx+1:end,:)*omegaks{1}(j)./sqrt(C(sub_indx+1:end,:)),S(sub_indx+1:end,:)*omegaks{2}(k).*sqrt(C(sub_indx+1:end,:)))./(C(sub_indx+1:end,:).*S(sub_indx+1:end,:));
            sup_filt{q} = trapz(scales{1}(sub_indx+1:end),trapz(scales{2},tmp,2));

    end
    sub_filt = cell2mat(reshape(sub_filt,length(omegaks{1}),length(omegaks{2})));
    sup_filt = cell2mat(reshape(sup_filt,length(omegaks{1}),length(omegaks{2})));
    
    %Calculate total energy
    total = sub_filt + sup_filt;
    C_delta = mean(real(total(:)));
    
    %Transform into physical domain
    sub_filt = ifftn(sub_filt); 
    sup_filt = ifftn(sup_filt);

    %Shift to center
    sub_filt = circshift(sub_filt,L(1:2)/2-1);
    sup_filt = circshift(sup_filt,L(1:2)/2-1);
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