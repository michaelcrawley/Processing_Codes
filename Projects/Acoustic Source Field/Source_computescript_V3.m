clear;
tic;
%Constants
downsample = 2;
FS = 4e5;
NCH = 18;
NFCH = 5:16;
D = 0.0254;
c = sqrt(1.4*287*(273.15+24));
rho = 1.2754;
observer.z = 100*D*cosd(30);
observer.r = 100*D*sind(30);
Nblocks = 1;
qmf = MakeONFilter('Battle',5);
multi = 2;

%Pick data and output directory
[piv_file,path] = uigetfile('*.mat','Identify NNE file');
load([path,filesep,piv_file]);
if ~exist('matversion','var'), matversion = 2; end
[acoustic_file,path] = uigetfile('*.bin','Identify Acoustic file');
fid = fopen([path,filesep,acoustic_file],'r');
outdir = 'St005_V5';
mkdir(outdir);

raw = fread(fid,'float32'); 
fclose(fid);
raw = reshape(raw,BS,NCH,[]);
raw = raw(:,:,1:Nblocks);
chk = y(1,:)>0;
r = y(:,chk)'/1000;
z = x(:,chk)'/1000;
[~,isort] = sort(r(:,1));
r = r(isort,:);
[M,N] = size(z);

%Find computation indices
indx = DS*width+1:downsample:BS-width*DS;
L = length(indx); 
t = (0:L-1)/(FS/downsample);
dt = mean(diff(t));
partial_t2 = mNumericalDerivative(2,2,1,L)/dt/dt;
field.c = c;
field.r = r;
field.z = z;
field.t = t;

tic;
for k = 1:Nblocks
    lafpa = raw(:,18,k);
    lafpa = lafpa(indx(1:L));
    inputs = raw(:,NFCH,k);
    disp(['Processing Block ',num2str(k),' of ',num2str(Nblocks)]);
    
    %Compute velocity field
    disp('Computing Velocity Field...');
    Uz = zeros(M,N,L);
    Ur = zeros(M,N,L);
    for n = 1:L
        tmp = inputs(indx(n)-DS*width:DS:indx(n)+DS*width,:);
        tmp = tmp(:)'/xt_norm;      
        switch matversion
            case 2  
                vel = nne(tmp)*d_norm;
                vel = reshape(vel,size(x,1),[]);
                Uz_tmp = vel(:,1:size(x,2))+ Um; 
                Ur_tmp = vel(:,size(x,2)+1:end) + Vm;
            case 3       
                vel = nne(tmp);
                vel = reshape(vel,size(x,1),[]);
                Uz_tmp = vel(:,1:size(x,2)).*d_norm + Um; 
                Ur_tmp = vel(:,size(x,2)+1:end).*d_norm + Vm;
            case 4
                ak = (nne(tmp).')*d_norm;
                vel = phi*ak;
                Uz_tmp = reshape(vel(1:numel(x)),size(x)) + Um;
                Ur_tmp = reshape(vel(numel(x)+1:end),size(x)) + Vm;    
        end
        Uz_tmp = Uz_tmp(:,chk)';
        Ur_tmp = Ur_tmp(:,chk)';
        Uz_tmp = Uz_tmp(isort,:);
        Ur_tmp = Ur_tmp(isort,:);
        Uz(:,:,n) = Uz_tmp;
        Ur(:,:,n) = Ur_tmp;  
    end
    clear phi;
    save([outdir,filesep,'Uz_blk',num2str(k),'.mat'],'Uz','lafpa','t','r','z');
    save([outdir,filesep,'Ur_blk',num2str(k),'.mat'],'Ur','lafpa','t','r','z');
    
    %Compute vortex ID
%     disp('Computing Vortex Identification...');
%     lambda = zeros(M,N,L);
%     parpool(8);
%     for n = 1:L
%         lambda(:,:,n) = SwirlingStrength(z,r,Uz(:,:,n),Ur(:,:,n));
%     end
%     delete(gcp);
%     save([outdir,filesep,'lambda_blk',num2str(k),'.mat'],'lambda','lafpa','t','r','z');
%     clear lambda;
    
    %Compute Solenoidal Velocity field
    disp('Computing Solenoidal Velocity Field...');
    s_Uz = zeros(M,N,L);
    s_Ur = zeros(M,N,L); 
    for n = 1:L
        [potential,solenoidal] = Helmholtz_Decomposition2DCyl_v4(z,r,Uz(:,:,n),Ur(:,:,n));
        s_Ur(:,:,n) = solenoidal.Ur;
        s_Uz(:,:,n) = solenoidal.Uz;
    end
    clear Uz Ur;
    save([outdir,filesep,'s_Uz_blk',num2str(k),'.mat'],'s_Uz','lafpa','t','r','z');
    save([outdir,filesep,'s_Ur_blk',num2str(k),'.mat'],'s_Ur','lafpa','t','r','z');
    
    %Compute PseudoPressure field
    disp('Computing Pseudo-Pressure Field...');
    ps = zeros(M,N,L);
    for n = 1:L
        ps(:,:,n) = Solenoidal_Pressure(z,r,rho,s_Uz(:,:,n),s_Ur(:,:,n));
    end
    clear s_Uz s_Ur;
    save([outdir,filesep,'ps_blk',num2str(k),'.mat'],'ps','lafpa','t','r','z');
    
    %Smooth in time
    disp('Smoothing...');
    sigma = std(ps(:));
    sps = reshape(ps,M*N,L);
    clear ps;
    T = zeros(1,2^nextpow2(L));
    parpool(8);
    parfor n = 1:M*N
        tmp = T;
        tmp(1:L) = sps(n,:);
        filt = ThreshWave(tmp,'S',1,sigma,multi,2,qmf);
        sps(n,:) = filt(1:L);
    end
    delete(gcp);
    sps = reshape(sps,[M,N,L]);
    save([outdir,filesep,'sps_blk',num2str(k),'.mat'],'sps','lafpa','t','r','z');
    
    %Compute Source
    disp('Computing Source Field...');
    source = zeros(M*N,L);
    sps = reshape(sps,M*N,L);
    for n = 1:M*N
        source(n,:) = partial_t2*sps(n,:)';
    end
    source = reshape(source,[M,N,L]);
    p = IntegrateRetardedTime(observer,field,source);
    save([outdir,filesep,'source_blk',num2str(k),'.mat'],'source','p','lafpa','t','r','z');
    clear source p sps;
end
toc