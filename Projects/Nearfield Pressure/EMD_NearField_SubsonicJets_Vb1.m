clc;
clear all;
close all;

CF = pwd; % Locate current folder

% for baseline jet
% PFN  = ([CF '\' 'Filtered Data' '\' 'BaselineJet']); 
% PFN  = ([CF '\' 'Filtered Data' '\' 'StDF = 0.02']); 
PFN = uigetdir('Specify Processing Directory:',pwd);

[FN] = getfiles('*.mat',PFN,'-s'); % Get the file name of processing data
ML = length(FN);

rD = [1.2 1.7 2.2 2.7 3.2 3.7 4.2 4.7 5.2 5.7 6.2 6.7 7.2 7.7 8.2];

for i = 1:1;
    if i == 1;
        IMF(i,:) = [ 2	2	2	2	2	2	2	3	3	3	3	4	4	4	4	4 ]; % rD 1.2, updated processing
%                                                         ||^||   
%     IMF(i,:) = [ 2  2  2	 2	 2	 2	 2   2	 3	 3	 3	 4	 4	 4	 4	 4  ]; % rD 1.2, first processing
    elseif i == 2;
        IMF(i,:) = [ 3	3	3	3	3	3	3	3	3	4	4	4	4	5	5	5 ]; % rD 1.7
    elseif i == 3;
        IMF(i,:) = [ 3	3	3	3	3	3	3	3	3	4	4	4	5	5	5	5 ];
    elseif i == 4;
        IMF(i,:) = [ 3	3	3	3	3	4	4	5	5	5	5	5	5	5	5	5 ];
    elseif i == 5;
        IMF(i,:) = [ 5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5 ];
    elseif i == 6;
        IMF(i,:) = [ 5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5 ];
    elseif i == 7;
        IMF(i,:) = [ 5	5	5	5	5	5	5	5	5	5	5	5	6	6	6	6 ];
    elseif i == 8;
        IMF(i,:) = [ 5	5	5	5	5	5	5	5	5	5	5	5	6	6	6	6 ];
    elseif i == 9;
        IMF(i,:) = [ 5	5	5	5	5	5	5	5	5	5	6	6	6	6	6	6 ];
    elseif i == 10;
        IMF(i,:) = [ 6	6	6	6	6	6	6	6	6	6	6	6	6	6	6	6 ];
    elseif i == 11;
        IMF(i,:) = [ 6	6	6	6	6	6	6	6	6	6	6	6	6	6	6	6 ];
    elseif i == 12;
        IMF(i,:) = [ 6	6	6	6	6	6	6	6	6	6	6	6	6	6	6	6 ];
    elseif i == 13;
        IMF(i,:) = [ 6	6	6	6	6	6	6	6	6	6	6	6	6	6	6	6 ];
    elseif i == 14;
        IMF(i,:) = [ 6	6	6	6	6	6	6	6	6	6	6	6	6	6	6	6 ];
    elseif i == 15;
        IMF(i,:) = [ 6	6	6	6	6	6	6	6	6	6	6	6	6	6	6	6 ];
    end;
end;

        
% IMF=[2 3 3 3 5 5 5 5 5 6 6 6 6 6 6]; % x/D = 2
% IMF=[2	3	3	3	5	5	5	5	5	6	6	6	6	6	6]; % x/D = 5
% IMF=[2	3	3	5	5	5	5	5	5	6	6	6	6	6	6]; % x/D = 8


% cpus = FindCoreCount
% matlabpool(cpus);

for si = 1:ML;
   
load([PFN filesep FN{ML-si+1}]); % load mat file from selected folder
FN{ML-si+1}  % FN{si}   

% NM = length(phys.x);
% [NSx,NSy,NSz]=size(nf.pblocks.p);
[NSx,NSy,NSz]=size(nf.pblocks.smp);

NM = phys.x; % selected mic position
LM = length(NM);

for i = 1:LM; % # of selected mic channels

%     tmp = nf.pblocks.p(:,i,:);
    tmp = nf.pblocks.smp(:,i,:);
    tmp = permute(tmp,[1 3 2]);
    p = reshape(tmp, [], 1);


size(p); % check input signal 

% Read in all experimental parameters
ao =phys.a;    % load ambient acoustic velocity
Uj =phys.Ue;   % load ambient acoustic velocity
FSt=phys.Stdf; % load forcing St 

% Parameter for experimental conditions
FS = 200000;
BS = 8192;
D  = 1*0.0254;
Uj = phys.Ue;
Md         = 1.0;
Mj         = phys.M;
gamma      = 1.4;

WNDO = 'hamming';

nfrq = ceil((BS+1)/2);  %# of distinct frequencies returned for SPL
f = (0:nfrq-1)'*(FS/BS);   %Frequency axis

if Mj < 1.0; 
    D_eff = D;
else; 
    D_eff = D*sqrt(Md/Mj)*((1+Mj*Mj*(gamma-1)/2)/(1+Md*Md*(gamma-1)/2))^((gamma+1)/(4*(gamma-1)));
end;

fc = Uj/D_eff;
St = f/fc;





    % Process based on the mic order 1 to 12
%     for mo = 6:6;
%             close all
            
            chtdata = strcat('ch1data = p(:,1);');
            eval(chtdata);
            
        % Signal conditioning
        m1 = mean(ch1data);
        ch1data = ch1data - m1;
%         m2 = mean(ch2data);
%         ch2data = ch2data - m2;
%         m3 = mean(ch3data);
%         ch3data = ch3data - m3;
%         m4 = mean(ch4data);
%         ch4data = ch4data - m4;
%         m5 = mean(ch5data);
%         ch5data = ch5data - m5;
%         m6 = mean(ch6data);
%         ch6data = ch6data - m6;
        
        yo = ch1data;
%         yo = ch1data';
        size(yo)
        
       
        % Read in the time series data
        MT = length (yo);
        dt = 0:1/FS:(MT-1)/FS;
        x  = dt';
        
        for RN = 1:7;

            if RN == 1;
                y = yo;
            else
                y  = yo - oy;
                yo = y;
            end;
            
            if RN == 7;
                yop = yo;
            end;

            % The iteration of sifting process for IMF
            DPN1 = 1;
            DPN2 = 1;
            fn   = 1
%             while (DPN1 > 0.05);
            while (DPN2 > 0.05);

                if DPN2 == 1;
                    y  = y;
                else
                    y  = y - my;
                    fn = fn + 1
                end;

                % plot the time series data
                figure(fn) 
                plot(x,y,'k-+')
                xlim([0.3,0.301]);


                % find positive peak 
                [py,fpx] = findpeaks(y);
                px = x(fpx);

                % plot positive peak envelope
                hold on;
                % plot(px,py,'bo');

                % construct an interpolating positive function using a cubic spline
                ipy = interp1(px,py,x,'spline');

                % check the end point for boundary condition
                if (ipy(1) < y(1));
                    ipy(1) = y(1);
                end;
                if (ipy(end) < y(end));
                    ipy(end) = y(end);
                end;

                plot(x,ipy,'b-');


                % find negative peak 
                [ny,fnx] = findpeaks(-1*y);
                nx = x(fnx);
                % reverse y vaules
                ny = -1*ny;

                % plot negative peak envelope
                hold on;
                % plot(nx,ny,'r-*');

                % construct an interpolating negative function using a cubic spline
                iny = interp1(nx,ny,x,'spline');

                % check the end point for boundary condition
                if (iny(1) > y(1));
                    iny(1) = y(1);
                end;
                if (iny(end) > y(end));
                    iny(end) = y(end);
                end;

                plot(x,iny,'r-');


                % Construct the average of the positive and negative function
                my = (ipy+iny)/2;
                size(my)
                plot(x,my,'g-','LineWidth',2);

                % Start the sifting process step 1
                SiP1 = 0.05;
                SiP2 = 10*SiP1;
                S1   = 0;
                S2   = 0;
                sig  = zeros(MT,1); 
                for kk = 1:MT;
                    sig1    = ipy(kk)+iny(kk);
                    sig2    = ipy(kk)-iny(kk);
                    if sig2 == 0;
                        sig(kk) = 0;
                    else
                        sig(kk) = sig1 / sig2;
                        % Sifting process step 2
                        if abs(sig(kk)) > SiP2;
                            S1 = 1 + S1;
                        end;
                        if abs(sig(kk)) > SiP1;
                            S2 = 1 + S2;
                        end;
                    end;
                end;
                DPN1 = S1/MT
                S1
                DPN2 = S2/MT
                S2

                plot(x,sig,'m--','LineWidth',2);
  
                clear px py nx ny fpx fnx ipy iny sig1 sig2 sig

            end; % while (DFN ==0);

%             % Document the convergence condition for DPN1 & DPN2
%             filename_DPN = strcat('ConvergenceCondition','_#',num2str(V_SN + Loop_Num -1),'a_Mic',num2str(mo),'.xlsx');
%             Out_DPN1 = strcat('C',num2str(RN+2));
%             xlswrite(filename_DPN,DPN1,'DPN',Out_DPN1);
%             Out_FN1 = strcat('D',num2str(RN+2));
%             xlswrite(filename_DPN,fn,'DPN',Out_FN1);
%             Out_DPN2 = strcat('E',num2str(RN+2));
%             xlswrite(filename_DPN,DPN2,'DPN',Out_DPN2);
%             Out_FN2 = strcat('F',num2str(RN+2));
%             xlswrite(filename_DPN,fn,'DPN',Out_FN2);
            
            
            oy = y;

%             cd IMF
%             % Output IMF to textfile
%             filename = strcat(char(V_FN),'_#',num2str(V_SN + Loop_Num -1),'a_Mic',num2str(mo),'_IMF_',num2str(RN),'.txt');
%             fid2 = fopen(filename,'wt');
%             fprintf(fid2,'%16.10f\t\n',y);
%             fclose(fid2);
%             cd ..
            
            
            % Output filtered raw time series data
            if RN == 7;
                OPY = strcat('ch',num2str(1+RN),'data = yop;');
                eval(OPY);
                
            else
                OPY = strcat('ch',num2str(1+RN),'data = y;');
                eval(OPY);
            end;

            yL = yo - y; % rest of the signal
            OPY = strcat('chLdata = yL;');
            eval(OPY);
            clear yL
            
            close all
            clc

        end; % for RN = 1:5;
    


% Check the addition effect w/ IMF4 + IMF5 + IMF6 + IMF7 for hydrodynamic field
ch9data = ch5data + ch6data + ch7data + ch8data;

% Check the addition effect w/ IMF1 + IMF2 + IMF3 for acoustic field
ch10data = ch2data + ch3data + ch4data;

if IMF(si,i) == 2;
    ch9data = ch4data + ch5data + ch6data + ch7data + ch8data;
    ch10data=ch2data + ch3data;
elseif IMF(si,i) == 3;
    ch9data = ch5data + ch6data + ch7data + ch8data;
    ch10data=ch2data + ch3data + ch4data;
elseif IMF(si,i) == 4;
    ch9data = ch6data + ch7data + ch8data;
    ch10data=ch2data + ch3data + ch4data + ch5data;
elseif IMF(si,i) == 5;
    ch9data = ch7data + ch8data;
    ch10data=ch2data + ch3data + ch4data + ch5data + ch6data;
elseif IMF(si,i) == 6;
    ch9data = ch8data;
    ch10data=ch2data + ch3data + ch4data + ch5data + ch6data + ch7data;
end;

hpms(si,i) = (rms(ch9data))^2;
apms(si,i) = (rms(ch10data))^2;
hp(si,i,:) = ch9data;
ap(si,i,:) = ch10data;

save 'StDF_0.02_Pms_hydro_aco_rD1.2.mat' hpms apms hp ap; % save data
% save 'baseline_Pms_hydro_aco.mat' hpms apms hp ap; % save data



for k = 1:10;
    ETP = strcat('p(:,k) = ch',num2str(k),'data;');
    eval(ETP);
end; 


% start to calculate auto-spectra
S = size(p);
if round(S(1)/BS)~=S(1)/BS
    error('Data is not an integer number of blocks')
end

MN = mean(p,1);
p = p - repmat(MN,S(1),1);   %subtracts mean from each channel - necessary to apply windowing
pM = reshape(p,BS,S(1)/BS*S(2));

%creates window as sparse diagonal matrix
if ~exist('WNDO','var')
    WNDO = 'rectwin';
end
WD = window(str2func(WNDO),BS);
WDM = spdiags(shiftdim(WD),0,BS,BS);
WDW = mean(WD.^2);  %Weight of windowing - must be scaled out of spectrum

rFFT = fft(WDM*pM,BS,1);  %Calculates FFT with the specified window
rFFT = rFFT(1:nfrq,:); %Discards symmetric portion of spectrum
rFFT = reshape(rFFT,nfrq,S(1)/BS,S(2));

mFFT = abs(rFFT).^2/(BS*FS*WDW);  %Scales spectrum according to sampling and weighting characteristics 
if rem(BS,2) %Conserves energy which would otherwise be lost by discardiing symmetric portion
    mFFT(2:end,:,:) = mFFT(2:end,:,:)*2;
else
    mFFT(2:end-1,:,:) = mFFT(2:end-1,:,:)*2;
end
SPL = squeeze(mean(mFFT,2)); %Calculates squared average of spectrum
size(SPL);
    
% SPL = 10.0*log10(fc)+ 10*log10(SPL/(0.00002^2));
SPL = SPL/(0.00002^2);
% SPL = SPL';

% BaSPL(si,:,:)=SPL; 

% end;
size(SPL)

figure(i)
semilogx(St,SPL(:,1:10),'LineWidth',2)
hold on;
set(gca,'FontWeight','bold','FontSize',12); % Set font as bold
xlabel('Strouhal Number','FontWeight','bold');
% ylabel('SPL per unit {\it{St}} (dB//(20{\it{\mu}}Pa^2))','FontWeight','bold');
ylabel('PSD','FontWeight','bold');
set(gca,'yscale','log');
xlim([0.01 6.5])
ylim([10^5 10^11])
set(gca,'xtick',[0.01 0.1 1 10]);
set(gca,'XTickLabel',{'0.01';'0.1';'1';'10'}); % Set the Xtick label
set(0,'DefaultLineLineWidth',2);
% set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');

% SFN = strcat('Baseline_EMD_Hydro_Acou_Spectra_rD',num2str(rD(si)),'_xD',num2str(phys.x(i),2));
% saveas(figure(i),[CF '\' 'EMD_Results' '\' 'Baseline' '\' SFN '.fig']);

SFN = strcat('StDF_0.15_EMD_Hydro_Acou_Spectra_rD',num2str(rD(si)),'_xD',num2str(phys.x(i),2));
saveas(figure(i),[CF '\' 'EMD_Results' '\' 'St 0.02' '\' SFN '.fig']);

% end; % for mo = 1:8;
end; %for i = 1:LM; % # of selected mic channels


end; % for si = 1:ML;
