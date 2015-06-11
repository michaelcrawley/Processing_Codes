% title('M = 1.3, T_o/T_a = 2.5, m = 0, St_{DF} = 0.6, CorrDist = 0.92\pm0.10 x/D')
% saveas(gcf,'WholeImageCorrelation-M13_TTR25_m0_St060.fig')
% saveFigure_v2(gcf,'WholeImageCorrelation-M13_TTR25_m0_St060',600)

%%%% SCHLIEREN %%%%%%%%%%%
IM = imrotate(imread('F:\Research\Mach1.3\Schlieren\20100928\AverageImages\vm11_St0.6_T475_ph000.png'),-2.75);
IM(:,:,2) = imrotate(imread('F:\Research\Mach1.3\Schlieren\20100928\AverageImages\vm11_St0.6_T475_ph090.png'),-2.75);
IM(:,:,3) = imrotate(imread('F:\Research\Mach1.3\Schlieren\20100928\AverageImages\vm11_St0.6_T475_ph180.png'),-2.75);
IM(:,:,4) = imrotate(imread('F:\Research\Mach1.3\Schlieren\20100928\AverageImages\vm11_St0.6_T475_ph270.png'),-2.75);
IM = double(IM)-128;
IM = IM(850:950,180:800,:);

dx = 1/160;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% PIV %%%%%%%%%%%%%%%%%
% IM = vm11_St099.CondAvg.Q; dx = mean(diff(m0_St028.x(:,1)));
% xl = find(m0_St028.x(:,1) > 8,1,'first');
% yl = find(m0_St028.y(1,:) < 0.75,1,'first');
% yl(2) = find(m0_St028.y(1,:) < -0.75,1,'first');
% IM = IM(1:xl,yl(1):yl(2),:); clear yl;
% 
% N = 3;
% IMt = interp2(IM(:,:,1),N); dx = dx*size(IM,1)/size(IMt,1);
% IMt(:,:,2) = interp2(IM(:,:,2),N);
% IMt(:,:,3) = interp2(IM(:,:,3),N);
% IMt(:,:,4) = interp2(IM(:,:,4),N); IM = IMt; clear IMt;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WHOLE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIM = 2;
L = size(IM,DIM);

clear corr
corr(1,:) = xcorr2_1d(IM(:,:,1),IM(:,:,1),DIM);
corr(2,:) = xcorr2_1d(IM(:,:,2),IM(:,:,1),DIM);
corr(3,:) = xcorr2_1d(IM(:,:,3),IM(:,:,1),DIM);
corr(4,:) = xcorr2_1d(IM(:,:,4),IM(:,:,1),DIM);

corr(5,:) = xcorr2_1d(IM(:,:,2),IM(:,:,2),DIM);
corr(6,:) = xcorr2_1d(IM(:,:,3),IM(:,:,2),DIM);
corr(7,:) = xcorr2_1d(IM(:,:,4),IM(:,:,2),DIM);

corr(8,:) = xcorr2_1d(IM(:,:,3),IM(:,:,3),DIM);
corr(9,:) = xcorr2_1d(IM(:,:,4),IM(:,:,3),DIM);

corr(10,:) = xcorr2_1d(IM(:,:,4),IM(:,:,4),DIM); clear IM;

Dx = (-(size(corr,2)-1)/2:(size(corr,2)-1)/2)*dx;

%% First Structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y0 = 476;   %pixels
% w = 100;
% B = 250;
% E = B+400;
% IMw = IM(y0-w:y0+w,B:E,:); clear IM;
% 
% clear corr
% corr(1,:) = xcorr2_1d(IMw(:,:,1),IMw(:,:,1),DIM);
% corr(2,:) = xcorr2_1d(IMw(:,:,2),IMw(:,:,1),DIM);
% corr(3,:) = xcorr2_1d(IMw(:,:,3),IMw(:,:,1),DIM);
% corr(4,:) = xcorr2_1d(IMw(:,:,4),IMw(:,:,1),DIM);
% 
% corr(5,:) = xcorr2_1d(IMw(:,:,2),IMw(:,:,2),DIM);
% corr(6,:) = xcorr2_1d(IMw(:,:,3),IMw(:,:,2),DIM);
% corr(7,:) = xcorr2_1d(IMw(:,:,4),IMw(:,:,2),DIM);
% 
% corr(8,:) = xcorr2_1d(IMw(:,:,3),IMw(:,:,3),DIM);
% corr(9,:) = xcorr2_1d(IMw(:,:,4),IMw(:,:,3),DIM);
% 
% corr(10,:) = xcorr2_1d(IMw(:,:,4),IMw(:,:,4),DIM);


%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(Dx,corr(1,:),'-b');
hold on
plot(Dx,corr(2,:),'-r');
plot(Dx,corr(3,:),'-g');
plot(Dx,corr(4,:),'-m');

plot(Dx,corr(5,:),'--b');
plot(Dx,corr(6,:),'--r');
plot(Dx,corr(7,:),'--g');

plot(Dx,corr(8,:),'-.b');
plot(Dx,corr(9,:),'-.r');

plot(Dx,corr(10,:),':b');

legend('0-0','0-90','0-180','0-270','90-90','90-180','90-270','180-180','180-270','270-270')
xlim([-4 4])
xlabel('Displacement (x/D)')
ylabel('Correlation Coefficient')
grid on

%% CALCULATE

% Auto-correlations
clear C
[M,I] = max(corr(1,round(L*1.1):end));
C(1) = I +round(L*0.1)-1;
[M,I] = max(corr(5,round(L*1.1):end));
C(5) = I +round(L*0.1)-1;
[M,I] = max(corr(8,round(L*1.1):end));
C(8) = I +round(L*0.1)-1;
[M,I] = max(corr(10,round(L*1.1):end));
C(10) = I +round(L*0.1)-1;

% Cross-correlations 1/4wavelength-forward
[M,I] = max(corr(2,L:end));
C(2) = (I-1)*4;
[M,I] = max(corr(6,L:end));
C(6) = (I-1)*4;
[M,I] = max(corr(9,L:end));
C(9) = (I-1)*4;

% Cross-correlations 1/2wavelength-forward
[M,I] = max(corr(3,L:end));
C(3) = (I-1)*2;
[M,I] = max(corr(7,L:end));
C(7) = (I-1)*2;

% Cross-correlations 3/4wavelength-forward
[M,I] = max(corr(4,L:end));
C(4) = round((I-1)*4/3);


% Cross-correlations 1/4wavelength-backward
[M,I] = max(corr(2,1:L));
C(11) = round((L-I)*4/3);
[M,I] = max(corr(6,1:L));
C(14) = round((L-I)*4/3);
[M,I] = max(corr(9,1:L));
C(16) = round((L-I)*4/3);

% Cross-correlations 1/2wavelength-backward
[M,I] = max(corr(3,1:L));
C(12) = (L-I)*2;
[M,I] = max(corr(7,1:L));
C(15) = (L-I)*2;

% Cross-correlations 3/4wavelength-backward
[M,I] = max(corr(4,1:L));
C(13) = (L-I)*4;

C = C*dx;
mean(C)
std(C)

% title('M = 1.3, T_o/T_a = 1.0, m = 0, St_{DF} = 0.26, CorrDist = 1.63\pm0.58 x/D')
% saveas(gcf,'WholeImageCorrelation-M13_TTR10_m0_St026.fig')
% saveFigure_v2(gcf,'WholeImageCorrelation-M13_TTR10_m0_St026',600)

