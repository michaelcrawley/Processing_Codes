function Wavelet_Decompose2DInst_V4(src_dir,flist,out_dir,optset,cpus)
    tic;
            
    %Set options
    if ~exist('optset','var'), optset = []; end
    if ~isfield(optset,'windowfun'), optset.windowfun = @rectwin; end
    if ~isfield(optset,'flag'), optset.flag = 'pblocks.smp'; end
    if ~isfield(optset,'mother'), optset.mother = 'morlet'; end
    if ~isfield(optset,'param'), optset.param = []; end
    if ~isfield(optset,'dj'), optset.dj = [1/24 1/24]; end
    if ~isfield(optset,'FS'), optset.FS = 2e5; end
    if ~isfield(optset,'Nx'), optset.Nx = 22; end
    if ~isfield(optset,'PAD'), optset.PAD = [0 3 0]; end
    if ~isfield(optset,'temp_flag'), optset.temp_flag = 'amb'; end
    if ~isfield(optset,'interp_method'), optset.interp_method = 'pchip'; end
    if ~isfield(optset,'D'), optset.D = 0.0254; end
    dt = 1/optset.FS;
    matversion = 1.2;
    fields = [1 2 3 7 14];
    
    %Open Matlab Pool if Necessary
    if ~exist('cpus','var')||isempty(cpus), cpus = 1; end
    if cpus > 1
        if matlabpool('size') > 0
            matlabpool close;
        end
        matlabpool(cpus);
    end
    
    for k = 1:length(flist)
        %load file
        tmpdata = load([src_dir filesep flist{k}]);
        [~,fname] = fileparts(flist{k});
        sig = eval(['tmpdata.nf.',optset.flag]);
        if strcmpi(optset.temp_flag,'amb')
            a = tmpdata.phys.a;
        elseif strcmpi(optset.temp_flag,'local')
            Temp = tmpdata.phys.To/(1+0.2*tmpdata.phys.M^2);
            a = sqrt(1.4*287*Temp);
        elseif strcmpi(optset.temp_flag,'jet')
            a = tmpdata.phys.Ue;
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
        
        if k == 1
            %Compute Analyzing Wavelets - scales{1} are temporal, scales{2}
            %are spatial
            [sub_filt,sup_filt,C_delta] = DaughterWavelets(L,[dt dx],optset.dj,a,optset.mother,optset.param);
            wndo_t = repmat(optset.windowfun(N(1)),[1 optset.Nx]); %window along t-dimension (DIM1)
        end
        
        %Compute Decompositions
        disp(['Processing File: ',fname,'...']);
        parfor n = 1:N(3)            
            %Interpolate 
            tmp = zeros(N(1),optset.Nx);
            for q = 1:N(1)
                tmp(q,:) = interp1(x,sig(q,:,n),int_x,optset.interp_method);
            end
            
            %%%%%%TEMPORARY
            %%%%%%TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tmp(:,13:end) = 0;
            
            
            intwaveform(:,:,n) = tmp;
            
            %Apply window, zero-pad            
            tmp = intwaveform(:,:,n).*wndo_t;
            swap = zeros(L(1),L(2));
            swap(1:N(1),1:optset.Nx) = tmp;
            tmp = swap;        
            
            %Filtered Reconstructions
            xh = fftn(tmp);
            tmp = ifftn(xh.*conj(sub_filt))/C_delta;
            subsonic(:,:,n) = real(tmp(:,1:optset.Nx));
            tmp = ifftn(xh.*conj(sup_filt))/C_delta;
            supersonic(:,:,n) = real(tmp(:,1:optset.Nx));
        end
        
        %Organize variables
        filter_params = optset;
        filter_params.a = a;
        filter_params.x = int_x;
        filter_params.file = mfilename;
        
        data.ff = ff;
        data.intwaveform = intwaveform;
        data.subsonic = subsonic;
        data.supersonic = supersonic;
        
        daq_params.phys = tmpdata.phys;
        daq_params.phys.D = optset.D;
        daq_params.filename = tmpdata.filename;
        daq_params.trigger = tmpdata.trigger;
        daq_params.phys.FS = optset.FS;
        
        %Get new filename
        [~,fname] = fileparts(flist{k});
        parts = regexp(fname,'_','split'); 
        short = parts(fields);
        newname =[sprintf('%s_',short{1:end-1}),short{end},' 2DWaveInst.mat'];
        
        %Save Output
        save([out_dir,filesep,newname],'data','filter_params','daq_params','matversion');
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
    if ~exist('param','var')||isempty(param), param = []; end
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
                    wavelet = @(omega,k) ((pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k > 0) + ... %positive wavenumbers
                                            (pi^-0.5)*exp(-0.5*((k + m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k < 0) ); %negative wavenumbers
                case 'morlet-improper'
                    if isempty(m), m = [6 6]; end
                    fcoef = m + sqrt(2 + m.^2);
                    wavelet = @(omega,k) ((pi^-0.5)*exp(-0.5*((k - m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k >= 0) + ... %positive wavenumbers
                                            (pi^-0.5)*exp(-0.5*((k + m(2)).^2)).*exp(-0.5*((omega - m(1)).^2)).*(k < 0) ); %negative wavenumbers
                case 'paul'
                    if isempty(m), m = 4; end
                    fcoef = (2*[m m]+1)/2; %need to find this number!!!
                    wavelet = @(omega,k) (2^m)/sqrt(m*factorial(2*m-1))*(((omega.^2+k.^2)/2).^m).*exp(-0.5*(omega.^2+k.^2)).*(omega > 0);
                otherwise
                    error('Undefined Mother Wavelet');
            end
        otherwise
            error('Incorrect Dimension Definition');
    end
end