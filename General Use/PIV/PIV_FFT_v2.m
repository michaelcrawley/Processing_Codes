function PIV_FFT_v2(x,y,U,V,Um,Vm,L1,L2,SW)

% L1 = 'Baseline - 155^o C';
% L2 = 'Baseline';
% SW = 1;

if SW
    Q = Um(round(length(y(1,:))/2),:);  %grab centerline profile
    CI = max(find( Q > max(Q)*0.9));    %find index for end of potential core
    CI = round(CI*0.85);    %reduces core length to ensure we are in core region
    clear Q
else
    CI = length(x(:,1));
end

FsX = abs(1/mean(diff(x(1:CI,1))));
FsY = abs(1/mean(diff(y(1,:))));
XR = max(x(1:CI,1))-min(x(1:CI,1));
YR = max(y(1,:))-min(y(1,:));

S = size(U);

mU = []; mV = []; pU = [];
for n = 1:S(3)
    progress(n,1,S(3),10)
    [kx,ky,mU(:,:,n),pU] = myfft2(FsX,FsY,U(1:CI,:,n)-Um(1:CI,:));
    [kx,ky,mV(:,:,n),pU] = myfft2(FsX,FsY,V(1:CI,:,n)-Vm(1:CI,:));
end
clear pU

mUavg = sqrt(mean(mU.^2,3));
mVavg = sqrt(mean(mV.^2,3));

figure
contourf(kx,ky,mUavg)
axis([1/XR max(kx(:)) 1/YR max(ky(:))]) %excludes contributions from wavelengths longer than viewing window
xlabel('k_{x/D} StreamWise')
ylabel('k_{y/D} Radial')
title({'RMS FFT of StreamWise (U) velocity fluctuations',['Mach0.9 - ' L1]})
colorbar
grid on
saveas(gcf,[cd L2 '_mU.fig']);
saveas(gcf,[cd L2 '_mU.png']);
figure
contourf(kx,ky,mVavg)
axis([1/XR max(kx(:)) 1/YR max(ky(:))]) %excludes contributions from wavelengths longer than viewing window
xlabel('k_{x/D} StreamWise')
ylabel('k_{y/D} Radial')
title({'RMS FFT of Radial (V) velocity fluctuations',['Mach0.9 - ' L1]})
colorbar
grid on
saveas(gcf,[cd L2 '_mV.fig']);
saveas(gcf,[cd L2 '_mV.png']);
figure
contourf(x(1:CI,:),y(1:CI,:),U(1:CI,:,1)-Um(1:CI,:))
xlabel('x/D StreamWise')
ylabel('y/D Radial')
title({'StreamWise (U) velocity fluctuations (m/s)',['Mach0.9 - ' L1]})
colorbar
grid on
saveas(gcf,[cd L2 '_u.fig']);
saveas(gcf,[cd L2 '_u.png']);
figure
contourf(x(1:CI,:),y(1:CI,:),V(1:CI,:,1)-Vm(1:CI,:))
xlabel('x/D StreamWise')
ylabel('y/D Radial')
title({'Radial (V) velocity fluctuations (m/s)',['Mach0.9 - ' L1]})
colorbar
grid on
saveas(gcf,[cd L2 '_v.fig']);
saveas(gcf,[cd L2 '_v.png']);

save([cd L2 '.mat'],'FsX','FsY','kx','ky','mUavg','mVavg','mU','mV')
