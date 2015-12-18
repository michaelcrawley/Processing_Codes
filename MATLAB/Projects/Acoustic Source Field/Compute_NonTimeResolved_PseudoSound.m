clear;

D = 0.0254;
rho = 1.2754;
src = '/mnt/Samimy_research/ACTIVE_DATA/Jet/Mach09/20150711/20150711/MAT';

%Find correct file
flist = getfiles('preprocessed.mat',src,'-a','-s');
St = flist{1}(5:9);
load([src filesep flist{1}]);

%Grab positive radii only and sort
r = data(1).Y/1000;
chk = r(1,:) >= 0;
r = r(:,chk);
z = data(1).X(:,chk)'/1000;
Uz = permute(data(1).U(:,chk,:),[2 1 3]);
Ur = permute(data(1).V(:,chk,:),[2 1 3]);
[~,Ir] = sort(r(1,:));
r = r(:,Ir)';
Uz = Uz(Ir,:,:);
Ur = Ur(Ir,:,:);

sol_Uz = zeros(size(Uz));
sol_Ur = sol_Uz;
sol_P = sol_Uz;
N = size(Uz);
for n = 1:N(3)
    [potential,solenoidal] = Helmholtz_Decomposition2DCyl_v2(z,r,Uz(:,:,n),Ur(:,:,n));
    sol_Uz(:,:,n) = solenoidal.Uz;
    sol_Ur(:,:,n) = solenoidal.Ur;
    sol_P(:,:,n) = Solenoidal_Pressure(z,r,rho,solenoidal.Uz,solenoidal.Ur);
end

save([St,'_Instantaneous Sol_UP.mat'],'-v7.3','src','flist','r','z','Uz','Ur','sol_Uz','sol_Ur','sol_P','badvec_chk');