function Wavelet_Decompose2DInst_V4(src_dir,flist,out_dir,a,flag,cpus)
%Computes STCWT in the physical domain, whereas V3 computed in the Fourier
%domain!

    tic;
    %Constants
    if ~exist('flag','var')||isempty(flag), flag = 'pblocks.p'; end
    if ~exist('cpus','var')||isempty(cpus), cpus = 1; end
    
    %Parameters
    mother = 'morlet';
    param = [6 6];
    dj = [1/12 1/12];
    dt = 1/2e5;
    D = 0.0254;
    Nx = 22;
    
    %Open Matlab Pool if Necessary
    if cpus > 1
        if matlabpool('size') ~= cpus
            if matlabpool('size') ~= 0
                matlabpool close;
            end
            matlabpool(cpus);
        end        
    end
    
    for k = 1:length(flist)
        %load file
        data = load([src_dir filesep flist{k}]);
        [~,fname] = fileparts(flist{k});
        sig = eval(['data.nf.',flag]);
                       
        %Get Parameters
        N = size(sig);
        t = (0:N(1)-1)*dt;
        x = data.phys.x*D; %convert to meters
        [X,T] = meshgrid(x,t); %for interpolating
        int_x = linspace(x(1),x(end),Nx); %interpolate onto regular grid
        dx = mean(diff(int_x));
        [new_X,new_T] = meshgrid(int_x,t);               
        
        
        %Grab any additional helpful data
        ff = eval(['data.ff.',flag]);
        
        %Initialize large matrices
        intwaveform = zeros(N(1),Nx,N(3));
        subsonic = intwaveform;
        supersonic = intwaveform;
        
        if k == 1
            %Compute Analyzing Wavelets - scales{1} are temporal, scales{2}
            %are spatial
            [sub_filt,sup_filt] = DaughterWavelets_physical([N(1),Nx,N(3)],[dt dx],dj,a,mother,param);
            filt_energy = abs(sub_filt).^2 + abs(sup_filt).^2;
        end
        
        %Compute Decompositions
        disp(['Processing File: ',fname,'...']);
        for n = 1:N(3)            
            %Interpolate
            intwaveform(:,:,n) = interp2(X,T,sig(:,:,n),new_X,new_T,'spline');
            
            qmax = N(1)*Nx;
            tmp_full = intwaveform(:,:,n);
            tmp_sub = zeros(N(1)*Nx,1);
            tmp_sup = tmp_sub;
            parfor q = 1:qmax
                [j,k] = ind2sub( [N(1), Nx], q );
                
                %Grab maps, compute energy
                sub_map = sub_filt((N(1)+1:2*N(1))-j,(Nx+1:2*Nx)-k);
                sup_map = sup_filt((N(1)+1:2*N(1))-j,(Nx+1:2*Nx)-k);
                c_delta = trapz(t,trapz(int_x,filt_energy((N(1)+1:2*N(1))-j,(Nx+1:2*Nx)-k),2));
                
                %subsonic                
                tmp = tmp_full.*sub_map;
                tmp_sub(q) = trapz(t,trapz(int_x,tmp,2))/c_delta;
                
                %supersonic
                tmp = tmp_full.*sup_map
                tmp_sup(q) = trapz(t,trapz(int_x,tmp,2))/c_delta;
            end
            
            %Filtered Reconstructions 
            subsonic(:,:,n) = reshape(tmp_sub,N(1),Nx);
            supersonic(:,:,n) = reshape(tmp_sup,N(1),Nx);
        end
        
        %Save Output
        save([out_dir,filesep,fname,' 2DWaveDecompPhys.mat'],'intwaveform','subsonic','supersonic','a','int_x','data','ff');
    end
    
    %Closes Parallel Pool, if Necessary
    if cpus > 1
        matlabpool('close');
    end
    disp(['Complete...',' Processing time: ',num2str(round(toc)),' seconds...']);
end

function [sub_filt,sup_filt] = DaughterWavelets_physical(L,dt,dj,ca,mother,param)
    %Grab Mother Wavelet info
    if ~exist('dj','var')||isempty(dj), dj = [0.125 0.125]; end
    if ~exist('mother','var')||isempty(mother), mother = 'morlet'; end
    if ~exist('param','var')||isempty(param), param = [6 6]; end
    if ~exist('ca','var')||isempty(ca), flag = false; else flag = true; end
    [~,fcoef] = MotherWavelets_physical('ST2',mother,param);
    
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
    
    %compute axes    
    t = (-(L(1)-1):(L(1)-1))*dt(1);
    x = (-(L(2)-1):(L(2)-1))*dt(2);
    L_ext = [length(t),length(x)];
    [C,S] = ndgrid(scales{1},scales{2}); 
    
    %Compute daughter wavelets
    wavelet = MotherWavelets_physical('ST2',mother,param);
    
    %Calculate subsonic-filtered map
    sub_filt = cell(L_ext);
    sup_filt = sub_filt;
    disp('Computing Phase Maps....');

    Nmax = prod(L_ext);
    parfor q = 1:Nmax            
            [j,k] = ind2sub( L_ext, q );
        
            %Subsonic map
            tmp = wavelet((t(j))./S(1:sub_indx,:).*sqrt(C(1:sub_indx,:)),(x(k))./S(1:sub_indx,:)./sqrt(C(1:sub_indx,:)))./(C(1:sub_indx,:).*S(1:sub_indx,:).^3);
            sub_filt{q} = trapz(scales{1}(1:sub_indx),trapz(scales{2},tmp,2));
            
            %Supersonic map
            tmp = wavelet((t(j))./S(sub_indx+1:end,:).*sqrt(C(sub_indx+1:end,:)),(x(k))./S(sub_indx+1:end,:)./sqrt(C(sub_indx+1:end,:)))./(C(sub_indx+1:end,:).*S(sub_indx+1:end,:).^3);
            sup_filt{q} = trapz(scales{1}(sub_indx+1:end),trapz(scales{2},tmp,2));
    end
    sub_filt = cell2mat(reshape(sub_filt,L_ext));
    sup_filt = cell2mat(reshape(sup_filt,L_ext));
end


function [wavelet,fcoef] = MotherWavelets_physical(dim,mother,m)
%This function returns a function handle, 'wavelet', given a dimension and
%string for the name of the mother wavelet. A parameter to modify the base
%mother wavelet can also be provided. For standard mother wavelets, the
%dimension provided will be an integer, for Spatio-Temporal mother
%wavelets, the dimension will be a string beginning with 'ST' and
%including the total number of dimensions.

    if ~exist('m','var'), m = []; end

    switch dim
        case 1
            error('Undefined Mother Wavelet');
        case 'ST2'
            switch lower(mother)
                case 'morlet'
                    if isempty(m), m = [6 6]; end
                    wavelet = @(t,x) (pi^-0.5)*exp(1i*(m(1)*x + m(2)*t)).*exp(-0.5*(x.^2 + t.^2));
                    fcoef = m + sqrt(2 + m.^2);
                otherwise
                    error('Undefined Mother Wavelet');
            end
        otherwise
            error('Incorrect Dimension Definition');
    end
end