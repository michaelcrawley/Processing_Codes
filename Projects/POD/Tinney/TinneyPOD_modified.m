%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author:         C. Tinney
%   Data Source:    Poitiers, France, (CoJEN)
%
%   Date Created:   March 8, 2006
%   Last Modififed: November 13, 2011
%
%   This will perform the 1-d POD on the line array data using a frequency 
%   dependent kernel.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%11111111111111111111111111111111111111111111111111111111111111111111111111
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%LineArray%%%%%LineArray%%%%%LineArray%%%%%LineArray%%%%%LineArray%%
%This will perform the 1-d POD on the line array using a frequencey
%dependent kernel.
% clear all; clc
% dirMfile = pwd;
% dirdat = '/users/cetinney/TINNEY/UT-Austin/UniversityInstruction/ASE396/Homework/HWPOD/PODData/';
% cd(dirdat);
font = 16;

X1 = 22;%20;                %Number of channels on the line array
N = 8192;%2^12;               %Sample size
B = 100;%20;                 %Total number of Blocks per position
fs = 2e5;%24096;             %Sampling frequency.
df = fs/N;              %frequency resolution
dt = 1/fs;              %time increment of data
T = N*dt;               %total time in each block of data [s]
De = 0.0254;%0.05;              %nozzle diameter [m]
Ue = 287;%125;               %exit velocity [m/s]
flab = [0:N-1]*(fs/N);  %frequency label [f]
tlab = [0:N-1]/fs*1000; %time label [ms]
StD = flab*De/Ue;       %Strouhal number
pref = 20e-6;           %standard reference pressure for air [Pa]
dx = 0.0254*cosd(8.6);%0.015;             %This is the intermicrophone spacing on the line array [m]
xlab = [dx:dx:X1*dx]/De;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This will create a kernel for the pressure line array: R(x,x,f)
Rppf(1:X1,1:X1,1:N) = 0;
Pvar = zeros(1,X1);
for b = 1:B, b
%     eval(strcat('load Pdatb',num2str(b+1000),'.mat'));
    pdata = data.intwaveform(:,:,b);
    pf = fft(pdata(1:N,1:X1));
    for j = 1:X1,
        for k = 1:X1,
            %This is the double-sided cross spectral density in units^2 and
            %is properly on account of Matlab's fft routine
            Rppf(j,k,:) = Rppf(j,k,:) + (T/N^2)*shiftdim(((conj(pf(:,j)).*pf(:,k))/B),-2);
        end
    end
    %This is the variance of the signal and does not use (N-1) in the
    %denominator...I have my reasons!...has to do with the ensemble
    %averaging process
    Pvar = Pvar + mean(pdata(:,1:X1).^2)/B;
end
% eval(strcat('save Rppf.mat Rppf'));
%return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lets look at some features of the data along the array
%This will filter out the auto power spectral density along the line array
for f = 1:N,
    Spp(f,:) = diag(Rppf(:,:,f))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the standard of deviation
Pstd = sqrt(Pvar)
%This is the integrated spectrum
Senergy(1:X1) = sum(Spp(:,:),1)*df
%This will determine the ratio between the integrated spectrum and the
%variance of the signal to check that it is appropriately scaled...uses
%Parseval's theorem. If the data is properly scaled it should have a ratio
%of 1.0!
ratio = (Pvar./Senergy)'
%This will correct the data if it is not properly scaled
for f = 1:N,
    Spp(f,:) = (Spp(f,:)).*ratio';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This will give you the autocorrelation function along the line array;
for x = 1:X1,
    rho(:,x) = (1/N)*real(ifft((N^2/T)*Spp(:,x)))/Pvar(x);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is in noise scaling (non-dimensional).  The factor two is because
%we want to view the single sided spectra.
for x = 1:X1,
    %This is the single sided spectra as a fraction of energy and without
    %1/f dependence.
    Gpp(:,x) = Spp(1:N/2,x).*flab(1:N/2)'/Pvar(x)*2;
    %This is the true fraction of energy but still has 1/f dependence
    Gpp2(:,x) = Spp(1:N/2,x)*df/Pvar(x)*2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since we are using units squared, we will use the 10*log10() formula.  If
%we took the square root of Saanew, then we would use the 20*log10()
%formula: i.e.Saadf(:,1:N) = 20*log10(sqrt(Spp')./pref);
SppdB(1:N,:) = (10*log10((Spp')./(pref^2)))';

%plot the PSD
font = 12; set(0,'DefaultLineMarkerSize',6.0); set(0,'DefaultLineLineWidth',1.5);
figure(1); loglog(flab(2:N/2),Spp(2:N/2,1),'k-'); hold on
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
loglog(flab(2:N/2),Spp(2:N/2,7),'r-',flab(2:N/2),Spp(2:N/2,14),'k--');
loglog(flab(2:N/2),Spp(2:N/2,20),'b--')
legend('mic 1','mic 7','mic 14','mic 20')
xlabel('f [Hz]'); ylabel('S_{pp}(f) [Pa^2/Hz]'); axis tight

figure(2); semilogx(StD(2:N/2),Gpp(2:N/2,1),'k-'); hold on
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
loglog(StD(2:N/2),Gpp(2:N/2,7),'r-',StD(2:N/2),Gpp(2:N/2,14),'k-');
loglog(StD(2:N/2),Gpp(2:N/2,20),'b--')
xlabel('St_D(f)'); ylabel('S_{pp}\cdotf/\sigma_p^2');
axis([StD(2) StD(N/2) 0 1.05])

figure(3); plot(tlab(1:N/8),rho(1:N/8,1),'k-'); hold on
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
loglog(tlab(1:N/8),rho(1:N/8,7),'r-',tlab(1:N/8),rho(1:N/8,14),'k-');
loglog(tlab(1:N/8),rho(1:N/8,20),'b--')
xlabel('St_D(f)'); ylabel('G_{pp}\cdotf/\sigma_p^2');
axis([tlab(1) tlab(N/8) -0.5 1.05]); grid on
%return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%this is the decomposition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,n,m] = size(Rppf);
phi(1:X1,1:X1,1:N) = 0;
for f = 1:N,
    %This will force the kernel's symmetry...you shouldn't have to do this,
    %but if you have negative Eigenvalues, this is usually why...has to do
    %with the acuracy of Matlab.  Also, a weighting function to replace the
    %numerical integration has not been applied.  This is ok to do so long
    %as the grid spacing is uniform...dx is the sme between measurement
    %points on the array.
    test = Rppf(:,:,f)';
    test = (test + Rppf(:,:,f))./2;
    Rppf(:,:,f) = test;
    
    %This will solve for the eigenvalues (Lambda) and eigenfunctions (Phi)
    [Phi,Lambda] = eig(Rppf(:,:,f));        %Solving the Eigenfunctions.
    lambda(1:n,f) = flipud(diag(Lambda));   %Flipping and taking the diagonal matrix.
    phi(:,:,f) = fliplr(Phi);
    
    Ltr(f) = trace(Lambda);
    energy_norm(1:n,f) = lambda(1:n,f)/Ltr(f);
end
Ltrsum = sum(Ltr);
energy = (lambda/Ltrsum)*100;

%This will check the scaling of the eigenvalues
Energy = sum(Spp,2)*dx;
figure(4); loglog(flab(2:N/2),Energy(2:N/2),'k',flab(2:N/2),Ltr(2:N/2)*dx,'r:')
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
legend('Integrated energy from PSD','Sum of Eigenvalues')
xlabel('f [Hz]'); ylabel('Energy')

%This will plot the eigenvalues
figure(5); semilogx(StD(1:N/2),(energy(1:4,1:N/2)*2)');
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
xlabel('St_{De}(f)'); ylabel('\Lambda^{(n)}(f)'); legend('n=1','n=2','n=3','n=4'); %grid on; %axis tight
axis([StD(1) StD(N/2) 0 5.0]); ylabel('\xi^{(n)}(f) [%]')
%return

%This will plot the eigenfunctions.  We will separate them into real and imaginary components.
for x = 1:X1,
    phireal(x,:) = real(phi(x,1,:));
    phiimag(x,:) = imag(phi(x,1,:));
end
[fx,fy] = ndgrid(flab,xlab);
figure(6); subplot(2,1,1); contourf(fx(5:200,:),fy(5:200,:),phireal(:,5:200)'); 
set(gca,'fontsize',font-2,'fontname','Times','fontangle','italic');
xlabel('f [Hz]'); ylabel('x/D'); axis([5 fx(200,1) xlab(1) xlab(X1)]); title('\phi(real)')
subplot(2,1,2); contourf(fx(5:200,:),fy(5:200,:),phiimag(:,5:200)');
set(gca,'fontsize',font-2,'fontname','Times','fontangle','italic');
xlabel('f [Hz]'); ylabel('x/D'); axis([5 fx(200,1) xlab(1) xlab(X1)]); title('\phi(imag)')

%eval(strcat('save Lambda.mat lambda'));
%eval(strcat('save Phi.mat phi'));
%clear phi; clear Rppf;
%return


%This will perform a low-dimensional reconstruction of the kernel
Bijtemp = zeros(X1,N);
for n = 1:X1, n
    for f = 1:N,
        Bij(:,f) =  lambda(n,f).*phi(:,n,f).*conj(phi(:,n,f));
    end
    Mij(:,:,n) = sqrt(real(Bij).^2+imag(Bij).^2);
    Bijtemp = Bijtemp + Mij(:,:,n);
end
figure; subplot(2,1,1); contourf(fx(5:200,:),fy(5:200,:),Mij(:,5:200,1)');
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
xlabel('f [Hz]'); ylabel('x/D'); title('R_f^{(1)}')
subplot(2,1,2); contourf(fx(5:200,:),fy(5:200,:),Mij(:,5:200,2)');
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
xlabel('f [Hz]'); ylabel('x/D'); title('R_f^{(2)}')

%I put this on a log scale so that the topography of the kernel can be more
%easily read.  These should be identical....to within the accuracy of
%Matlab.
figure; subplot(2,1,1); contourf(fx(5:200,:),fy(5:200,:),log(Bijtemp(:,5:200,1))');
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
xlabel('f [Hz]'); ylabel('x/D'); title('Fully low-dimensional reconstruction')
colorbar
subplot(2,1,2); contourf(fx(5:200,:),fy(5:200,:),log(Spp(5:200,:)));
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
xlabel('f [Hz]'); ylabel('x/D'); title('Original PSD')
colorbar
%return

%This will check to see that the eigenvalues are properly scaled

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This will create a kernel for the pressure line array: R(x,x,f)
af(1:X1,1:N) = 0;
anm(1:X1,1:N) = 0;
for b = 1:B, b
%     eval(strcat('load Pdatb',num2str(b+1000),'.mat'));
    pdata = data.intwaveform(:,:,b);
    pf = fft(pdata(1:N,1:X1))';
    for n = 1:X1,
        for f = 1:N,
            a(n,f) = sum(pf(:,f).*conj(phi(:,n,f)));
        end
        %This may be complex due to matlab accuracy errors.  Let's take the
        %real part just in case.
        anm(n,:) = anm(n,:) + real((T/N^2)*a(n,:).*conj(a(n,:))/B);
    end
end

figure(20);
loglog(flab(2:N/2),lambda(1,2:N/2),'k',flab(2:N/2),anm(1,2:N/2),'r:'); hold on
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
loglog(flab(2:N/2),lambda(3,2:N/2),'k-',flab(2:N/2),anm(3,2:N/2),'b:');
xlabel('f [Hz]'); ylabel('y');
legend('\lambda^{(1)}(f)','|{a^{(3)}(f)}^2|','\lambda^{(3)}(f)','|{a^{(3)}(f)}^2|')
axis tight



%%%%%%%%%%%%%%%%%%this is the reconstruction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This will perform the first projection using only one of the time series
pf = (fft(pdata))';
%This will calculate the random POD expansion coefficients
for n = 1:X1, n
    for f = 1:N,
        a(n,f) = sum(pf(:,f).*conj(phi(:,n,f)));
    end
end
%This will reconstruct a low dimensional image of the time series
%First we start with the POD modes (n), and they are in descending order
n1 = 1; n2 = 4;
for j = 1:X1,
    for f = 1:N,
        uPOD(j,f) = 0;
        uPOD1(j,f) = a(1,f)*phi(j,1,f);
        uPOD2(j,f) = a(2,f)*phi(j,2,f);
        uPOD3(j,f) = a(3,f)*phi(j,3,f);
        uPOD4(j,f) = a(4,f)*phi(j,4,f);
        for n = n1:n2,
            uPOD(j,f) = uPOD(j,f) + a(n,f)*phi(j,n,f);
        end
    end
end
pPOD = real(ifft(uPOD'));
pPOD1 = real(ifft(uPOD1'));
pPOD2 = real(ifft(uPOD2'));
pPOD3 = real(ifft(uPOD3'));
pPOD4 = real(ifft(uPOD4'));

Pstd = sqrt(mean(pdata.^2));
PstdPOD = sqrt(mean(pPOD.^2));
PstdPOD1 = sqrt(mean(pPOD1.^2));
PstdPOD2 = sqrt(mean(pPOD2.^2));
PstdPOD3 = sqrt(mean(pPOD3.^2));
PstdPOD4 = sqrt(mean(pPOD4.^2));

set(0,'DefaultLineMarkerSize',4.0); set(0,'DefaultLineLineWidth',1.5);
figure; plot(tlab(1:400),pdata(1:400,10),'k-o',tlab(1:400),pPOD(1:400,10),'b'); hold on
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
plot(tlab(1:400),pPOD1(1:400,10),'r',tlab(1:400),pPOD2(1:400,10),'k--');
plot(tlab(1:400),pPOD3(1:400,10),'k-.',tlab(1:400),pPOD4(1:400,10),'k.');
legend('raw signal','n = 1:4 (4/20)','n = 1','n = 2','n = 3','n = 4',1); xlabel('t [s]'); ylabel('p(t) [Pa]')
axis tight

set(0,'DefaultLineMarkerSize',6.0);
figure; plot(xlab,Pstd,'k',xlab,PstdPOD,'bo',xlab,PstdPOD1,'r^',xlab,PstdPOD2,'kp'); hold on
plot(xlab,PstdPOD3,'k<',xlab,PstdPOD4,'k*');
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
legend('raw signal','n = 1:4','n = 1','n = 2','n = 3','n = 4',2); xlabel('x/D'); ylabel('\sigma_p [Pa]')
axis tight


[Xloc Tloc]=meshgrid(xlab,tlab/1000*Ue/De);
figure; subplot(2,1,1); contourf(Tloc(1:1200,:),Xloc(1:1200,:),pdata(1:1200,:),10);
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
title('original time series'); xlabel('t*U_e/D_e'); ylabel('x/D'); axis tight
subplot(2,1,2); contourf(Tloc(1:1200,:),Xloc(1:1200,:),pPOD(1:1200,:),10);
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
title('Filtered low-dimensional reconstruction using n = 1:4'); xlabel('t*U_e/D_e'); ylabel('x/D'); axis tight

figure; subplot(2,1,1); contourf(Tloc(1:1200,:),Xloc(1:1200,:),pPOD1(1:1200,:),10);
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
title('Filtered low-dimensional reconstruction using n = 1'); xlabel('t*U_e/D_e'); ylabel('x/D'); axis tight
subplot(2,1,2); contourf(Tloc(1:1200,:),Xloc(1:1200,:),pPOD2(1:1200,:),10);
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
title('Filtered low-dimensional reconstruction using n = 2'); xlabel('t*U_e/D_e'); ylabel('x/D'); axis tight

figure; subplot(2,1,1); contourf(Tloc(1:1200,:),Xloc(1:1200,:),pPOD3(1:1200,:),10);
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
title('Filtered low-dimensional reconstruction using n = 3'); xlabel('t*U_e/D_e'); ylabel('x/D'); axis tight
subplot(2,1,2); contourf(Tloc(1:1200,:),Xloc(1:1200,:),pPOD4(1:1200,:),10);
set(gca,'fontsize',font,'fontname','Times','fontangle','italic');
title('Filtered low-dimensional reconstruction using n = 4'); xlabel('t*U_e/D_e'); ylabel('x/D'); axis tight

%Observations:  POD mode 1 appears to be capturing the evanescent
%hydrodynamic signatures that are convecting at subsonic speeds. Whereas
%POD mode 3 appears to be capturing acoustic (sonic) signatures attempting
%to propogate to the far-field regions.
