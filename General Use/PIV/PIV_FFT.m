clear mU mW pU pW

FsX = abs(1/mean(diff(x(:,1))));
FsZ = abs(1/mean(diff(z(1,:))));
XR = max(x(:,1))-min(x(:,1));
ZR = max(z(1,:))-min(z(1,:));

S = size(U);

mU = []; mW = []; pU = []; pW = [];
for n = 1:S(3)
    progress(n,1,S(3),10)
%     [kx,kz,mU(:,:,n),pU(:,:,n)] = myfft2(FsX,FsZ,U(:,:,n)-Um);
%     [kx,kz,mW(:,:,n),pW(:,:,n)] = myfft2(FsX,FsZ,W(:,:,n)-Wm);
    [kx,kz,mU(:,:,n),pU] = myfft2(FsX,FsZ,U(:,:,n)-Um);
    [kx,kz,mW(:,:,n),pW] = myfft2(FsX,FsZ,W(:,:,n)-Wm);
end

mUavg = sqrt(mean(mU.^2,3));
mWavg = sqrt(mean(mW.^2,3));
% pUavg = sqrt(mean(pU.^2,3));
% pWavg = sqrt(mean(pW.^2,3));

figure
contourf(kx,kz,mUavg)
axis([1/XR max(kx(:)) 1/ZR max(kz(:))]) %excludes contributions from wavelengths longer than viewing window
xlabel('k_{x/D} StreamWise')
ylabel('k_{z/D} CrossStream')
title({'RMS FFT of StreamWise (U) velocity fluctuations','Mach0.9 - Vm11 - 3.5kHz - 38^o C'})
colorbar
grid on
saveas(gcf,[cd '\1\Vm11_U.fig']);
figure
contourf(kx,kz,mWavg)
axis([1/XR max(kx(:)) 1/ZR max(kz(:))]) %excludes contributions from wavelengths longer than viewing window
xlabel('k_{x/D} StreamWise')
ylabel('k_{z/D} CrossStream')
title({'RMS FFT of CrossStream (W) velocity fluctuations','Mach0.9 - Vm11 - 3.5kHz - 38^o C'})
colorbar
grid on
saveas(gcf,[cd '\1\Vm11_W.fig']);
figure
contourf(x,z,U(:,:,1)-Um)
xlabel('x/D StreamWise')
ylabel('z/D CrossStream')
title({'StreamWise (U) velocity fluctuations (m/s)','Mach0.9 - Vm11 - 3.5kHz - 38^o C'})
colorbar
grid on
saveas(gcf,[cd '\1\Vm11_U_fluct.fig']);
figure
contourf(x,z,W(:,:,1)-Wm)
xlabel('x/D StreamWise')
ylabel('z/D CrossStream')
title({'CrossStream (W) velocity fluctuations (m/s)','Mach0.9 - Vm11 - 3.5kHz - 38^o C'})
colorbar
grid on
saveas(gcf,[cd '\1\Vm11_W_fluct.fig']);

clear CL TKE S n pU pW
save([cd '\1\Vm11.mat'])
