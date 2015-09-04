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
observer.z = 22*D;
observer.r = 9*D;

%Load data
load('FFNBP_arch64_St005_UV.mat')
fid = fopen('/mnt/Samimy_research/ACTIVE_DATA/Jet/Mach09/20150711/Acoustic/M09_St005_Sync_T23.4.bin','r');

% load('FFNBP_St025_arch64_UV.mat')
% fid = fopen('/mnt/Samimy_research/ACTIVE_DATA/Jet/Mach09/20150711/Acoustic/M09_St025_Sync_T24.1.bin','r');
raw = fread(fid,'float32'); 
fclose(fid);
raw = reshape(raw,BS,NCH,[]);
NF = raw(:,NFCH,1);
lafpa = raw(:,18,1);
clear raw;
chk = y>0;
r = reshape(y(chk),size(x,1),[])'/1000;
z = reshape(x(chk),size(x,1),[])'/1000;
r = flipud(r);
[M,N] = size(z);
clear xt;

%Find computation indices
indx = DS*width+1:downsample:BS-width*DS;
L = length(indx); 
t = (0:L-1)/(FS/downsample);
dt = mean(diff(t));
lafpa = lafpa(indx(1:L));

%Initialize variables
lighthill = zeros(M,N,L);
p_s = lighthill;
Uz = lighthill;
Ur = lighthill;
p_Uz = lighthill;
s_Uz = lighthill;

%Compute time-resolved velocity --> solenoidal velocity --> incomp source
parpool(12);
parfor n = 1:L
%     counter = counter + 1;
%     disp(['Processing interation: ',num2str(counter),'/',num2str(L)]);
    %grab/format inputs, compute estimated velocity
    tmp = inputs(indx(n)-DS*width:DS:indx(n)+DS*width,:);
    tmp = tmp(:)'/xt_norm;        
    vel = nne(tmp)*d_norm;
    vel = reshape(vel,size(x,1),[]);
    Uz_tmp = vel(:,1:size(x,2)) + Um; 
    Ur_tmp = vel(:,size(x,2)+1:end) + Vm;
        
    %Compute Solenoidal velocity field
    Uz_tmp = flipud(reshape(Uz_tmp(chk),N,[])');
    Ur_tmp = flipud(reshape(Ur_tmp(chk),N,[])');
    Uz(:,:,n) = Uz_tmp;
    Ur(:,:,n) = Ur_tmp;    
    [potential,solenoidal] = Helmholtz_Decomposition2DCyl_v2(z,r,Uz_tmp,Ur_tmp);
    p_Uz(:,:,n) = potential.Uz;
    s_Uz(:,:,n) = solenoidal.Uz;
    
    %Compute and filter source
    S = ComputeAASource(r,z,rho,solenoidal.Ur,solenoidal.Uz);
    S(1:2,:) = 0;
    S(end-1:end,:) = 0;
    S(:,1:2) = 0;
    S(:,end-1:end) = 0;
    Ssm = nanmoving_average2(S,3,3);
    lighthill(:,:,n) = Ssm;
    
    p_s(:,:,n) = Solenoidal_Pressure(z,r,rho,solenoidal.Uz,solenoidal.Ur);
end
% 
save('Uz.mat','Uz','z','r','t','lafpa');
save('Ur.mat','Ur','z','r','t','lafpa');
save('p_Uz.mat','p_Uz','z','r','t','lafpa');
save('s_Uz.mat','s_Uz','z','r','t','lafpa');
save('lighthill.mat','lighthill','z','r','t','lafpa');
save('p_s.mat','p_s','z','r','t','lafpa');

toc
% exit
% partial_t2 = mNumericalDerivative(2,2,1,L)/dt/dt;
% for m = 1:M
%     for n = 1:N
%         tmp = p_s(m,n,:);
%         source(m,n,:) = partial_t2*tmp(:);
%     end
% end
% 
source = zeros(M*N,L);
signal = reshape(smoothed,[M*N L]);
parfor k = 1:(M*N)
    tmp = signal(k,:);
    source(k,:) = partial_t2*tmp(:);
end
source = reshape(source,[M N L]);

% spatial_cutout = repmat(tukeywin(N,.05).',[M 1]).*exp(-(2e5)*(r.^4));

% 
% 
% 
% %Compute observer pressure
field.c = c;
field.r = r;
field.z = z;
field.t = t;
% p = IntegrateRetardedTime(observer,field,lighthill);