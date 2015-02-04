function SoundProcess_v8(M,AT,AP,RH)
%-v8 will accept directories containing baselines whose file name uses a
%forcing frequency of value zero. The original files with form
%(Mxx_xxxx_F0.0_xxxx_Txxxx_xxxx.NOS) will be renamed as
%(Mxx_Baseline_Txxxx.NOS) before the file list is presented to the user.
%-v7 modifies OASPL calculations to use non-smoothed spectra as much as 
%possible and to compute dOASPL for both detoned and non-detoned forced
%spectra if user requested detoning.
%-v6 adds the correction for atmospheric absorption. It updates the first 
%user call to verify hard-coded parameters. It corrects the method for
%calculating the FFT. Previous versions were calculating SPL levels which
%were dependent on sampling properties.
%-v5 adds the correction for microphone orientation.
%
%INPUTS
%  M - Mach number
%  AT - Ambient Temperature (C)
%  AP - Ambient Pressure (kPa)
%  RH - Relative Humidity (%)
%
%For the processing of acoustic data.  Program produces several files of
%processing results for an arbitrary number of channels of data.  Only one
%data directory may be processed per program execution.  Program requires
%initial (not automatable) user input to specify file locations and other
%parameters.  The program works on the assumption that there are N blocks
%of data in a file and that a RMS-FFT should be created from those N
%blocks.  It is assumed that any one directory being processed was
%collected at a constant ambient temperature.
%
%Code assumes isentropically/ideally expanded flow.
%
%The program cannot tell if files have been previously processed and will
%overwrite existing output files.
%
%    REQUIREMENTS
% - All files must have .NOS extension in order to be processed.
% - All files must be tab delimited and contain only data.
% - All data for a given set of parameters must be in a single file (i.e.
% no spreading channels among multiple files).
% - Data must be arranged so that each column represents a channel.
% - Data must have an integer number of blocks (i.e. if there are 1000
% samples per block and 10 blocks, there better be 10000 rows in the file -
% no more, no less).  One consequence of this is that all channels must
% have the same sampling rate, samples per block, and number of samples.
% - All data files must have "_Txxx" in the file name where xxx is a
% floating point value for the jet stagnation temperature in Celcius (xxx 
% can be more than three characters).
% - If forcing harmonics are to be removed, file must have "_Fxxx" in the
% file name where xxx is a floating point value for the forcing frequency 
% in kHz (xxx can be more than three characters).
% - Baseline and forcing cases must be in the same directory. A baseline
% case should have the following form: "Mxx_Baseline_Txxx.NOS".
% - File naming should follow this kind of convention
% "M0.9_m0_F010_T243.1_DAUH.NOS" - The important thing is that each parameter be
% separated by an underscore. This example was for: Mach 0.9, mode zero
% forcing, 10 kHz forcing frequency, 243.1 Celcius jet stagnation
% temperature, and automatic duty cycle forcing.
% - A Corresponding .NOS calibration file for each channel must exist. 
% These files must all be in one directory, but it doesn't have to be the 
% same directory as the data.  These files must contain one column of data
% or N columns in ascending channel order where N is the number of channels.
% - All calibration files must have "_CHx" in there name where x is the
% number of the channel (x can be more than one character).
% - Calibration files should be noise data acquired when microphone is
% exposed to tone generator signal.  Program will assume that the peak
% value in the calibration spectrum is this reference signal.  The
% parameter ref is the amplitude (dB) of this signal.
% - Microphone distance measurements (x/D) must be known.  These will be
% manually inputted by the user during execution.
% - Ambient temperature (C), pressure (kPa), and relative humidity (%) must
% be known.  It is assumed that single values will suffice for all data 
% being processed.
%
%    OUTPUT FILES
% - .fftNOS - File containing: Frequency axis, Strouhal number axis, and a
% column of SPL (dB) data for each channel acquired.
% - .HR.fftNOS - File containing: Frequency axis, Strouhal number axis, and
% a column of SPL (dB) data in which harmonics have been removed and profile 
% has been smoothed for each channel acquired.
% - .fig - MATLAB figure containing plots of SPL spectra for each channel
% acquired.
% - .png - standard picture file containing copy of image in .fig file.
% - .SET - Tab delimited file containing useful documented/calculated 
% parameters for the data.


    %Hard-coded parameters which may need to be changed
% AT = 28;    %Ambient Temperature - C
% AP = 100.5; %Atmospheric Pressure - kPa
% RH = 50;    %Relative Humidity - %
% M = 1.3;   %Mach Number
D = 0.0254; %jet diameter - m
ref = 114;  %reference amplitude - dB
sampleRate = 200000;    %data sampling rate - Hz
dSize = 8192;   %Samples per block
NormD = 80; %Normalization distance - x/D
ChPol = [90 80 70 60 50 45 40 35 30 25];   %Mic polar angles - degrees
MO = [0 10 20 30 40 45 50 55 60 0];     %Mic orientation angle relative to incident sound - degrees - Note: MO==0 is normal incidence
ChD = [49 49.75 52.25 56.5 64 69.25 76.25 85.5 98 116]; %Mic distance - r/D
DT = 1;     %0 = No tone subtraction (DeToning)
dBShift = 0;   %db Offset of spectra - this arbitrary offset can be used to correct spectra for uncalibrated amplification levels
DPoints = 6;   %Trim first and last N points off spectra.
WNDO = 'hamming';    %Data windowing method - use "rectwin" for no window
SD = 'F:\Research\';   %Starting directory
	
   
   %Asks user to verify proper parameters
disp('Current Processing Values: ');
disp(['                Ambient Temp (C): ' num2str(AT)]);
disp(['          Ambient Pressure (kPa): ' num2str(AP)]);
disp(['           Relative Humidity (%): ' num2str(RH)]);
disp(['                Jet Diameter (m): ' num2str(D)]);
disp(['                     Mach Number: ' num2str(M)]);
disp(['     Reference Signal Level (dB): ' num2str(ref)]);
disp(['              Sampling Rate (Hz): ' num2str(sampleRate)]);
disp(['               Samples Per Block: ' num2str(dSize)]);
disp(['Mic Normalization Distance (r/D): ' num2str(NormD)]);
disp(['      Mic Polar Angles (degrees): ' num2str(ChPol,'%.0f ')]);
disp([' Mic Orientation Angle (degrees): ' num2str(MO,'%.0f ')]);
disp(['              Mic Distance (r/D): ' num2str(ChD,'%.2f ')]);
if DT==0
    disp('                Tone Subtraction: NO');
else
    disp('                Tone Subtraction: YES');
end
disp(['   Spectrum Vertical Offset (dB): ' num2str(dBShift)]);
disp(['      Trim ends of sprectrum (N): ' num2str(DPoints)]);
disp(['                     Data Window: ' num2str(WNDO)]);
disp(['              Starting Directory: ' SD]);
yn = input('  Is this information correct? (y/n): ','s');
if ~strcmpi(yn,'y')
    error('Fix Values in code and run again')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FILE SELECTION/BASIC CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Asks user for a directory containing data - can only process one
    %directory per program execution.
dir_name = uigetdir(SD,'Specify Directory Containing Data'); 
flist = struct2cell(dir(dir_name));
q = strfind(flist(1,:),'.NOS'); %extracts list of data files - ignores all other files
keep = ~cellfun('isempty',q);
flist = flist(1,logical(keep));
flist = IDBaselines(dir_name,flist);    %rewrites file names of Baselines
[keep,ok] = listdlg('PromptString','Select files to process:',...
                'SelectionMode','multiple',...
                'ListString',flist);    %Asks the user to select the subset of data files they wish to process
if ok==0
    error('Program Terminated Due to User Selection of Cancel')
end
flist = flist(keep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALIBRATION FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Asks user to select the calibration files to use and verifies proper
    %naming convention.
[FileName,PathName] = uigetfile('*.NOS','Select Calibration File(s) for Data',SD,'MultiSelect','on');

%%%%% All operations beyond this point are done without user interaction. %%%%%

Nch = length(FileName); %Determines the number of channels of data by the number of calibration files
q = strfind(FileName, '_CH');
if isempty(q{1})
    error('Channel Number Prefix Must be All Caps...CH not Ch')
end
ChN = zeros(1,Nch);
for n = 1:Nch   %Extracts channel number from from file name
    qq = FileName{n}(q{n}+3:q{n}+8);
    qqq = strfind(qq,'_');
    ChN(n) = str2num(qq(1:qqq(1)-1));
end
[B,IX] = sort(ChN); clear B;
FileName = FileName(IX); clear IX; %Reorders files in ascending channel number to ensure proper correlation to data columns

if length(ChPol) ~= Nch
    error('    Microphone information mismatch!! Number of calibration files does not match microphone location information.')
end

for n = 1:Nch   %Calculates spectral reference value for each channel from the calibration files
    tmp = dlmread([PathName FileName{n}],'\t');
    if size(tmp,2)==1   %If there is only one channel of data, assume it is for the nth channel
        rRaw(:,n) = tmp;
    else    %Else, assume that the nth column contains the calibration for the nth channel
        rRaw(:,n) = tmp(:,n);
    end
end
[fx,SPL] = calcSPL_v2(rRaw,sampleRate,dSize,WNDO);
fx = fx(DPoints:end-DPoints); %First and last points will be thrown out as garbage
ChRef = sum(SPL(2:end-1,:),1)*sampleRate/dSize;   %calculates OASPL of calibration signal. This OASPL == prms^2
ChRef = spdiags(1./ChRef(:),0,Nch,Nch); %prepares the calibration values for scaling the data.

ChDM = spdiags(20*log10(ChD(:)/NormD) +dBShift,0,Nch,Nch); %prepares normalization distance for scaling the data.

color = 'brgkymc';  %Initializes array of colors for plotting

    %Inserts extra \ into directory path so it will print properly
Q = fliplr(strfind(PathName,'\'));
for nn = 1:length(Q)
    PathName = [PathName(1:Q(nn)) PathName(Q(nn):end)];
end

disp('  Calibration Data Processed...Beginning Data Processing')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:length(flist)
    Fnm = flist{n};
    disp(['    Processing File: ' Fnm])  %Displays file currently being processed
        %Locates temperature information in file name
    TLoc = strfind(Fnm,'_T');   
    Tend = min(strfind(Fnm(TLoc+1:end),'_'))+TLoc;
    if isempty(Tend)
        Tend = length(Fnm)-length('.NOS')+1;
    end
    Ttxt = Fnm(TLoc+2:Tend-1);
    Ttxt(strfind(Ttxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
    To = str2num(Ttxt)+273.15;  %Converts to Kelvin
        %Locates forcing frequency information in file name if it exists
    FLoc = strfind(Fnm,'_F');   
    if ~isempty(FLoc)
        SHarmonics=true;    %This boolean will be used to determine if forcing frequency information needs to be written in output file
        Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
        FH = str2num(Fnm(FLoc+2:Fend-1))*1000; clear FLoc Fend;
    else
        SHarmonics=false;   
        FH = 0;
    end
    
    Temp = To/(1+1/5*M^2);  %Calculates isentropically expanded jet exit temperature - K
    Std = fx*D/(M*sqrt(7/5*287.05*Temp));    %Creates Strouhal number axis assuming Mach number is constant
    ReN = 1.01e5/(287.05*Temp) *M*sqrt(7/5*287.05*Temp)*...
        D /(1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2));  %Calculates jet exit Reynolds number using Sutherlands formula for viscosity and ideal gas law for density
    
    rRaw = dlmread([dir_name filesep Fnm],'\t'); %reads data file
    [f,SPL] = calcSPL(rRaw,sampleRate,dSize,WNDO); clear rRaw f;
    
    dB = ref+ 10*log10(SPL*ChRef);    %Converts into dB using calibration and reference values
    dB = dB(DPoints:end-DPoints,:);   %removes junk points at beginning and end
    
    dBCal = dB+ ones(size(dB))*ChDM;     %Normalizes data to a standard distance
    dBCal = micOrientationCorrect(fx,MO,dBCal); %Corrects spectra for microphone orientation 
        %Corrects for absorption by atmosphere - Spectra are now representative 
        %of sound propagated to microphone distance with no absorption.
    alpha = atmAbsorption(fx,AT,AP,RH);
    [x,y] = meshgrid(ChD*D,alpha);
    dBCal = dBCal + x.*y;
    clear alpha x y
    
    dBCal2 = dBCal;
    if DT~=0
        for m = 1:Nch   %Smoothes and removes forcing harmonics from spectrum. If FH==0 (forcing frequency) program doesn't look for harmonics
            dBCal2(:,m) = ToneSubtract_v4(fx,dBCal(:,m),FH,10);   %Old function call
%             dBCal2(:,m) = ToneSubtract_v5(fx,dBCal(:,m),'Poly',30);
        end
    else
        for m = 1:Nch   %Smoothes forcing harmonics from spectrum.
            dBCal2(:,m) = ToneSubtract_v4(fx,dBCal(:,m),0,10);   %Old function call
%             dBCal2(:,m) = ToneSubtract_v5(fx,dBCal(:,m),'Poly',30);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OUTPUT FILE GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%  PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure  %Plots all spectra on a single figure
    semilogx(Std,dBCal(:,1),color(2))
    LBL = 'Ch 1';  %Variable contains labels which will become column headers on output files
    hold on
    for m = 2:Nch
        semilogx(Std,dBCal(:,m),color(mod(m,7)+1))
        LBL = [LBL '\tCh ' num2str(m)];
    end
    LEGN = mat2cell(ChPol,1,ones(1,Nch));
    for m = 1:Nch   %Plots smoothed spectra on top of original spectra
        semilogx(Std,dBCal2(:,m),color(mod(m-1,7)+1))
        LEGN{m} = num2str(LEGN{m});
    end
    legend(LEGN,'Location','NorthWest')
    xlabel(['Strouhal Number for D = ' num2str(D) ' m Jet at M = ' num2str(M)])
    ylabel('SPL (dB)')
    title({['Average Spectrum for: ' Fnm],['Data Normalized for ' num2str(NormD) ' x/D'],['Reynolds Number: ' num2str(ReN,'%.0f')]})
    grid on
    xlim([0.01 4])
    saveFigure_v2(gcf,[dir_name filesep Fnm(1:end-4)],300)   %Saves figure as .png
    saveas(gcf,[dir_name filesep Fnm(1:end-4) '.fig'])   %Saves figure as .fig
    close    %Closes figure
    
    %%%%%  SPECTRAL DATA FILES  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([dir_name filesep Fnm(1:end-3) 'fftNOS'],'w');   %Opens .fftNOS file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number\t' LBL '\n']);
    fclose(fid);
    dlmwrite([dir_name filesep Fnm(1:end-3) 'fftNOS'],[fx Std dBCal], 'delimiter', '\t', '-append');   %writes spectral data to file
    
    fid = fopen([dir_name filesep Fnm(1:end-3) 'S.fftNOS'],'w');   %Opens .HR.fftNOS file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number\t' LBL '\n']);
    fclose(fid);
    dlmwrite([dir_name filesep Fnm(1:end-3) 'S.fftNOS'],[fx Std dBCal2], 'delimiter', '\t', '-append');   %writes smoothed spectral data to file

    %%%%%  SET FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([dir_name filesep Fnm(1:end-3) 'SET'],'w');
    fprintf(fid,'PROCESS VALUES FOR DATA\n');
    fprintf(fid,['Sampling Rate (Hz):\t' num2str(sampleRate) '\n']);
    fprintf(fid,['Block Size:\t' num2str(dSize) '\n']);
    fprintf(fid,['Reference Amplitude (dB):\t' num2str(ref) '\n']);
    fprintf(fid,['Window:\t' WNDO '\n']);
    fprintf(fid,['Normalization Distance (x/D):\t' num2str(NormD) '\n']);    
    fprintf(fid,['Number of Channels:\t' num2str(Nch) '\n']);
        %Records location and name of each calibration file used.
    fprintf(fid,['Location of Calibration Files:\t' PathName '\n']);
    fprintf(fid,'Calibration Files:\t');
    for nn = 1:length(FileName)
        fprintf(fid,[FileName{nn} '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mic Distance (x/D):\t');
    for nn = 1:length(ChD)
        fprintf(fid,[num2str(ChD(nn)) '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mic Polar Angle (deg):\t');
    for nn = 1:length(ChPol)
        fprintf(fid,[num2str(ChPol(nn)) '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mic Orientation Angle (deg):\t');
    for nn = 1:length(MO)
        fprintf(fid,[num2str(MO(nn)) '\t']);
    end
    fprintf(fid,'\n');
    if DT==0
        fprintf(fid,'Tone Subtraction:\tNO\n');
    else
        fprintf(fid,'Tone Subtraction:\tYES\n');
    end
    fprintf(fid,['Spectrum offset (dB):\t' num2str(dBShift) '\n']);
    fprintf(fid,['Trim Ends of Spectrum (N):\t' num2str(DPoints) '\n']);
    fprintf(fid,['Jet Diameter (m):\t' num2str(D) '\n']);
    fprintf(fid,['Mach Number:\t' num2str(M) '\n']);
    fprintf(fid,['Ambient Temperature (K):\t' num2str(AT+273.15) '\n']);
    fprintf(fid,['Ambient Pressure (kPa):\t' num2str(AP) '\n']);
    fprintf(fid,['Relative Humidity (%%):\t' num2str(RH) '\n']);
    fprintf(fid,['Stagnation Temperature (K):\t' num2str(To) '\n']);
    fprintf(fid,['Exit Temperature (K):\t' num2str(Temp) '\n']);
    fprintf(fid,['Total Temperature Ratio:\t' num2str(To/(AT+273.15)) '\n']);
    fprintf(fid,['Exit Velocity (m/s):\t' num2str(M*sqrt(7/5*287.05*Temp)) '\n']);
    fprintf(fid,['Jet Velocity Ratio:\t' num2str(M*sqrt(Temp/(AT+273.15))) '\n']);
    fprintf(fid,['Exit Density (kg/m^3):\t' num2str(1.01e5/(287.05*Temp)) '\n']);
    fprintf(fid,['Exit Viscosity (Pa s):\t' num2str((1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2))) '\n']);   %equation is Sutherlands formula - 291.15 is reference temperature
    fprintf(fid,['Jet Reynolds Number:\t' num2str(ReN) '\n']);
    
        %Calculates and prints OASPL for all channels
    fprintf(fid,'OASPL (dB):\t');    
    Bg = min(find(Std > 0.04)); %Excludes spectral components below 0.04 Strouhal number
    Ed = min(find(Std > 4));    %Excludes spectral components above 4 Strouhal number
    OASPL = 10*log10(trapz(fx(Bg:Ed),10.^(dBCal(Bg:Ed,:)/10),1)); clear dBCal;
    for nn = 1:length(ChD)
        fprintf(fid,[num2str(OASPL(nn)) '\t']);
    end
    fprintf(fid,'\n');
        %Calculates and prints the average acoustic energy
    fprintf(fid,'Average Energy (dB):\t');
    delta_Pol = abs(diff(ChPol)); IW = ([delta_Pol/2 0] + [0 delta_Pol/2])/(max(ChPol)-min(ChPol));
    AAE = 10*log10(sum(IW.*10.^(OASPL/10)));
    fprintf(fid,[num2str(AAE) '\n']);
        %Calculates and prints the peak SPL magnitude and Strouhal number
    fprintf(fid,'Peak of SPL Curve\n');    
    fprintf(fid,'  Magnitude (dB):\t');
    for nn = 1:length(ChD)
        [c,I] = max(dBCal2(1:Ed,nn));
        I(I < 6) = 6;
        c = mean(dBCal2(I-5:I+5,nn));
        fprintf(fid,[num2str(c) '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'  Frequency (Std):\t');
    for nn = 1:length(ChD)
        [c,I] = max(dBCal2(1:Ed,nn));
        fprintf(fid,[num2str(Std(I)) '\t']);
    end
    fprintf(fid,'\n');
        %Calculates and prints the high frequency slope of the spectra 
        %including uncertainty as seen on semilogx scale.
    fprintf(fid,'Slope of SPL Curve (Strouhal axis):\t');
    for nn = 1:length(ChD)
        [c,I] = max(dBCal2(1:Ed,nn));   %Ignores all points to left of peak
        Pval = mmpolyfit(log(Std(I*2:Ed)),dBCal2(I*2:Ed,nn),1,'Weight',1./Std(I*2:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
            %Calculates uncertainty of fit
        Serr = sum((Pval(1)*log(Std(I*2:Ed))+Pval(2)-dBCal2(I*2:Ed,nn)).^2);
        merr = length(Std(I*2:Ed));
        Derr = merr*sum(log(Std(I*2:Ed)).^2) - sum(log(Std(I*2:Ed)))^2;
        PE = sqrt(Serr*merr/(merr-2)/Derr);
        
        fprintf(fid,[num2str(Pval(1),'%10.3f') ' +/- ' num2str(PE,'%10.3f') '\t']);
    end
    fprintf(fid,'\n');
    
        %If forcing is present, calculates/prints several relevant values
    if SHarmonics
        fprintf(fid,['Forcing Frequency (Hz):\t' num2str(FH) '\n']);
        fprintf(fid,['Forcing Strouhal Number:\t' num2str(FH*D/(M*sqrt(7/5*287.05*Temp))) '\n']);
            %Calculates and prints OASPL for all channels
        if DT~=0
            fprintf(fid,'OASPL_Detoned (dB):\t');
            OASPL_DT = 10*log10(trapz(fx(Bg:Ed),10.^(dBCal2(Bg:Ed,:)/10),1));
            for nn = 1:length(ChD)
                fprintf(fid,[num2str(OASPL_DT(nn)) '\t']);
            end
            fprintf(fid,'\n');
                %Calculates and prints the average acoustic energy
            fprintf(fid,'Average Energy_Detoned (dB):\t');
            AAE_DT = 10*log10(sum(IW.*10.^(OASPL_DT/10)));
            fprintf(fid,[num2str(AAE_DT) '\n']);
        end
        
            %Finds SET file(s) for Baseline.  If none are found, the
            %following calculations are not performed.
        slist = struct2cell(dir(dir_name));
        qs = strfind(slist(1,:),'Baseline');
        qs2 = strfind(slist(1,:),'.SET');
        kp = zeros(size(qs));
        for m = 1:length(qs)
            if (~isempty(qs{m}))&&(~isempty(qs2{m}))
                kp(m) = 1;
            end
        end
        if sum(kp) > 0  %If baseline SET file(s) are found
            slist = slist(1,logical(kp));
            OASPLB = zeros(1,Nch);
            AAEB = 0;
            OASPLB_DT = OASPLB;
            AAEB_DT = 0;
            for nn = 1:length(slist)    %reads baseline SET file(s) and extracts OASPL and total measured energy
                [NOASPLB,VB,kp] = readSETfile([dir_name filesep slist{nn}],[28 29 33 34]);
                OASPLB(nn,:) = VB{1};
                AAEB(nn) = VB{2};
                OASPLB_DT(nn,:) = VB{3};
                AAEB_DT(nn) = VB{4};
            end
            OASPLB = mean(OASPLB,1);    %averages OASPL from multiple baselines
            AAEB = mean(AAEB);  %averages average acoustic energy from multiple baselines
            OASPLB_DT = mean(OASPLB_DT,1);    %averages OASPL from multiple baselines
            AAEB_DT = mean(AAEB_DT);  %averages average acoustic energy from multiple baselines
            
                %Calculates and prints delta OASPL
            fprintf(fid,'dOASPL (dB):\t');      
            for nn = 1:length(ChD)
                fprintf(fid,[num2str(OASPL(nn) - OASPLB(nn)) '\t']);
            end
            fprintf(fid,'\n');
                %Calculates and prints change in average acoustic energy
            fprintf(fid,'dAAE (dB):\t');      
            fprintf(fid,[num2str(AAE-AAEB) '\n']);
            
            if DT~=0
                    %Calculates and prints delta OASPL
                fprintf(fid,'dOASPL_Detoned (dB):\t');      
                for nn = 1:length(ChD)
                    fprintf(fid,[num2str(OASPL_DT(nn) - OASPLB_DT(nn)) '\t']);
                end
                fprintf(fid,'\n');
                    %Calculates and prints change in average acoustic energy
                fprintf(fid,'dAAE_Detoned (dB):\t');      
                fprintf(fid,[num2str(AAE_DT-AAEB_DT) '\n']);
            end
        end       
    else
            %Calculates and prints OASPL for all channels
        fprintf(fid,'OASPL_Smoothed (dB):\t');    
        OASPL_DT = 10*log10(trapz(fx(Bg:Ed),10.^(dBCal2(Bg:Ed,:)/10),1)); clear dBCal;
        for nn = 1:length(ChD)
            fprintf(fid,[num2str(OASPL_DT(nn)) '\t']);
        end
        fprintf(fid,'\n');
            %Calculates and prints the average acoustic energy
        fprintf(fid,'Average Energy_Smoothed (dB):\t');
        AAE_DT = 10*log10(sum(IW.*10.^(OASPL_DT/10)));
        fprintf(fid,[num2str(AAE_DT) '\n']);
    end
    fclose(fid);
end
disp(' Done')   %Displays "Done" when all files have been processed

