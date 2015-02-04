function SoundProcess_v9(M,AT,AP,RH)
%-v9
%  - Saves all results into a directory whose name is a date-time stamp 
%created at run-time.
%  - Produces a function dependency report and creates an archive 
%containing: all m-files used to process the data, the calibration files, 
%and a list of which files were processed. 
%  - Baseline data is now linearly interpolated versus temperature to more 
%accurately match the characteristics at any given forced case.
%  - Smoothed spectra are now computed as a 30th order polynomial fit which
%works better than moving average detoning for narrowband peaks (see
%ToneSubtract_v5.m).
%  - Computes the Normalized Jet Power (NJP) and writes that additional
%information to SET files.
%  - Computes correction for atmospheric with switch for choosing lossless
%propagation or standard day propagation.
%  - Computes normalization distance propagation with switch for
%propagation correction which accounts for a noise source distribution.
%See function "NormDistance.m".
%  - Dynamically scans baseline SET files for locations of needed
%information instead of assuming known locations.
%  - Properly ignores additional filename parts (e.g. Mxx_EXAMPLE_xxxx.NOS)
%so that additional file identifiers may be used.
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
NDs = 'Radial';     %Normalization distance correction method - inputs: 'Radial'-simple radial correction; 'NS'-noise source location correction
NormD = 80; %Normalization distance - x/D
ChPol = [90 80 70 60 50 45 40 35 30 25];   %Mic polar angles - degrees
MO = [0 10 20 30 40 45 50 55 60 0];     %Mic orientation angle relative to incident sound - degrees - Note: MO==0 is normal incidence
ChD = [49 49.75 52.25 56.5 64 69.25 76.25 85.5 98 116]; %Mic distance - r/D
DT = 1;     %0 = No tone subtraction (DeToning)
dBShift = 0;   %db Offset of spectra - this arbitrary offset can be used to correct spectra for uncalibrated amplification levels
DPoints = 6;   %Trim first and last N points off spectra.
WNDO = 'hamming';    %Data windowing method - use "rectwin" for no window
atmAbs = 'Standard';    %Atmospheric absorption correction - inputs: 'Standard'-correct to standard day; 'Lossless'-removes absorption
SD = 'Z:\My Documents\Projects\';   %Starting directory

pref = 20e-6;   %Pascals - reference pressure - used in calculating normalized jet power
   
   %Asks user to verify proper parameters
disp('Current Processing Values: ');
disp(['                 Ambient Temp (C): ' num2str(AT)]);
disp(['           Ambient Pressure (kPa): ' num2str(AP)]);
disp(['            Relative Humidity (%): ' num2str(RH)]);
disp(['                 Jet Diameter (m): ' num2str(D)]);
disp(['                      Mach Number: ' num2str(M)]);
disp(['      Reference Signal Level (dB): ' num2str(ref)]);
disp(['               Sampling Rate (Hz): ' num2str(sampleRate)]);
disp(['                Samples Per Block: ' num2str(dSize)]);
disp(['Mic Normalization Distance Method: ' NDs]);
disp([' Mic Normalization Distance (r/D): ' num2str(NormD)]);
disp(['       Mic Polar Angles (degrees): ' num2str(ChPol,'%.0f ')]);
disp(['  Mic Orientation Angle (degrees): ' num2str(MO,'%.0f ')]);
disp(['               Mic Distance (r/D): ' num2str(ChD,'%.2f ')]);
if DT==0
    disp('                 Tone Subtraction: NO');
else 
    disp('                 Tone Subtraction: YES');
end
disp(['    Spectrum Vertical Offset (dB): ' num2str(dBShift)]);
disp(['       Trim ends of sprectrum (N): ' num2str(DPoints)]);
disp(['                      Data Window: ' num2str(WNDO)]);
disp(['                   Atm Absorption: ' atmAbs]);
disp(['               Starting Directory: ' SD]);
yn = input('  Is this information correct? (y/n): ','s');
if ~strcmpi(yn,'y')
    error('Fix Values in code and run again')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FILE SELECTION/BASIC CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Asks user for a directory containing data - can only process one
    %directory per program execution.
src_dir = uigetdir(SD,'Specify Directory Containing Data'); 
flist = struct2cell(dir(src_dir));
q = strfind(flist(1,:),'.NOS'); %extracts list of data files - ignores all other files
keep = ~cellfun('isempty',q);
flist = flist(1,logical(keep));
flist = IDBaselines(src_dir,flist);    %rewrites file names of Baselines
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
[CalFiles,CalPath] = uigetfile('*.NOS','Select Calibration File(s) for Data',SD,'MultiSelect','on');

%%%%% All operations beyond this point are done without user interaction. %%%%%

Nch = length(CalFiles); %Determines the number of channels of data by the number of calibration files
q = strfind(CalFiles, '_CH');
if isempty(q{1})
    error('Channel Number Prefix Must be All Caps...CH not Ch')
end
ChN = zeros(1,Nch);
for n = 1:Nch   %Extracts channel number from from file name
    qq = CalFiles{n}(q{n}+3:q{n}+8);
    qqq = strfind(qq,'_');
    ChN(n) = str2double(qq(1:qqq(1)-1));
end
[B,IX] = sort(ChN); clear B;
CalFiles = CalFiles(IX); clear IX; %Reorders files in ascending channel number to ensure proper correlation to data columns

if length(ChPol) ~= Nch
    error('    Microphone information mismatch!! Number of calibration files does not match microphone location information.')
end

for n = 1:Nch   %Calculates spectral reference value for each channel from the calibration files
    tmp = dlmread([CalPath CalFiles{n}],'\t');
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


color = 'brgkymc';  %Initializes array of colors for plotting

    %Inserts extra \ into directory path so it will print properly
Q = fliplr(strfind(CalPath,'\')); CalPathP = CalPath;
for nn = 1:length(Q)
    CalPathP = [CalPathP(1:Q(nn)) CalPathP(Q(nn):end)];
end

disp('  Calibration Data Processed...Beginning Data Processing')

tstamp = now;
out_dir = [src_dir filesep datestr(tstamp,30)];
mkdir(out_dir);


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
    To = str2num(Ttxt)+273.15;  %#ok<ST2NM> %Converts to Kelvin
        %Locates forcing frequency information in file name if it exists
    FLoc = strfind(Fnm,'_F');   
    if ~isempty(FLoc)
        SHarmonics=true;    %This boolean will be used to determine if forcing frequency information needs to be written in output file
        Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
        FH = str2num(Fnm(FLoc+2:Fend-1))*1000; clear FLoc Fend; %#ok<ST2NM>
    else
        SHarmonics=false;   
        FH = 0;
    end
    
    Temp = To/(1+1/5*M^2);  %Calculates isentropically expanded jet exit temperature - K
    Std = fx*D/(M*sqrt(7/5*287.05*Temp));    %Creates Strouhal number axis assuming Mach number is constant
    ReN = 1.01e5/(287.05*Temp) *M*sqrt(7/5*287.05*Temp)*...
        D /(1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2));  %Calculates jet exit Reynolds number using Sutherlands formula for viscosity and ideal gas law for density
    
    rRaw = dlmread([src_dir filesep Fnm],'\t'); %reads data file
    [f,SPL] = calcSPL_v2(rRaw,sampleRate,dSize,WNDO); clear rRaw f;
    
    dB = ref+ 10*log10(SPL*ChRef);    %Converts into dB using calibration and reference values
    dB = dB(DPoints:end-DPoints,:);   %removes junk points at beginning and end
    
    dBCal = micOrientationCorrect(fx,MO,dB); %Corrects spectra for microphone orientation 
        %Corrects for absorption by atmosphere - Spectra are now representative 
        %of sound propagated to microphone distance with no absorption.
    [alpha,alpha_s] = atmAbsorption_v2(fx,AT,AP,RH,atmAbs); 
    [x,y] = meshgrid(ChD*D,alpha); dBCal = dBCal + x.*y;    %Removes absorption at current conditions
    clear alpha x y
    
    if strcmpi('Radial',NDs)    %Correct to normalized distance using simple radial propagation
        [dBCal, al] = NormDistance(dBCal,ChD,NormD,dBShift);
        if strcmpi('Standard',atmAbs)   %Adds absorption at standard day conditions if selected
            [x,y] = meshgrid(ChD*D.*al',alpha_s); dBCal = dBCal - x.*y;
        end
    else        %Correct to normalized distance using noise source location
        [dBCal, al] = NormDistance(dBCal,Std,ChD,ChPol*pi/180,NormD,dBShift);
        if strcmpi('Standard',atmAbs)
            [x,y] = meshgrid(ChD*D,alpha_s); dBCal = dBCal - x.*y.*al;
        end
    end
    clear al x y alpha_s
    
    dBCal2 = dBCal;
    if DT~=0
        for m = 1:Nch   %Smoothes and removes forcing harmonics from spectrum. If FH==0 (forcing frequency) program doesn't look for harmonics
            dBCal2(:,m) = ToneSubtract_v5(fx,dBCal(:,m),'Poly',30);
        %If tones are broadband and poly fitting is failing, use this alternate call:
%             dBCal2(:,m) = ToneSubtract_v5(fx,dBCal(:,m),'Moving',FH,10);
        end
    else
        for m = 1:Nch   %Smoothes forcing harmonics from spectrum.
            dBCal2(:,m) = ToneSubtract_v5(fx,dBCal(:,m),'Poly',30);
        %If tones are broadband and poly fitting is failing, use this alternate call:
%             dBCal2(:,m) = ToneSubtract_v5(fx,dBCal(:,m),'Moving',0,10);
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
        LBL = [LBL '\tCh ' num2str(m)]; %#ok<AGROW>
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
    saveas(gcf,[out_dir filesep Fnm(1:end-4) '.fig'])   %Saves figure as .fig
    saveFigure_v2(gcf,[out_dir filesep Fnm(1:end-4)],300)   %Saves figure as .png
    close    %Closes figure
    
    %%%%%  SPECTRAL DATA FILES  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([out_dir filesep Fnm(1:end-3) 'fftNOS'],'w');   %Opens .fftNOS file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number\t' LBL '\n']);
    fclose(fid);
    dlmwrite([out_dir filesep Fnm(1:end-3) 'fftNOS'],[fx Std dBCal], 'delimiter', '\t', '-append');   %writes spectral data to file
    
    fid = fopen([out_dir filesep Fnm(1:end-3) 'S.fftNOS'],'w');   %Opens .HR.fftNOS file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number\t' LBL '\n']);
    fclose(fid);
    dlmwrite([out_dir filesep Fnm(1:end-3) 'S.fftNOS'],[fx Std dBCal2], 'delimiter', '\t', '-append');   %writes smoothed spectral data to file

    %%%%%  SET FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([out_dir filesep Fnm(1:end-3) 'SET'],'w');
    fprintf(fid,'PROCESS VALUES FOR DATA\n');
    fprintf(fid,['Processing Time Stamp:\t' datestr(now,0) '\n']);
    fprintf(fid,['Sampling Rate (Hz):\t' num2str(sampleRate) '\n']);
    fprintf(fid,['Block Size:\t' num2str(dSize) '\n']);
    fprintf(fid,['Reference Amplitude (dB):\t' num2str(ref) '\n']);
    fprintf(fid,['Window:\t' WNDO '\n']);
    fprintf(fid,['Normalization Distance (x/D):\t' num2str(NormD) '\n']);    
    fprintf(fid,['Normalization Method:\t' NDs '\n']);
    fprintf(fid,['Number of Channels:\t' num2str(Nch) '\n']);
        %Records location and name of each calibration file used.
    fprintf(fid,['Location of Calibration Files:\t' CalPathP '\n']);
    fprintf(fid,'Calibration Files:');
    for nn = 1:Nch
        fprintf(fid,['\t' CalFiles{nn}]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mic Distance (r/D):'); fprintf(fid,'\t%.2f',ChD); fprintf(fid,'\n');
    fprintf(fid,'Mic Polar Angle (deg):'); fprintf(fid,'\t%.2f',ChPol); fprintf(fid,'\n');
    fprintf(fid,'Mic Orientation Angle (deg):'); fprintf(fid,'\t%.2f',MO); fprintf(fid,'\n');
    if DT==0
        fprintf(fid,'Tone Subtraction:\tNO\n');
    else
        fprintf(fid,'Tone Subtraction:\tYES\n');
    end
    fprintf(fid,['Spectrum offset (dB):\t' num2str(dBShift) '\n']);
    fprintf(fid,['Trim Ends of Spectrum (N):\t' num2str(DPoints) '\n']);
    if strcmpi('Standard',atmAbs)
        fprintf(fid,'Atmospheric Absorption:\tStandard Day\n');
    else
        fprintf(fid,'Atmospheric Absorption:\tLossless Propagation\n');
    end
    fprintf(fid,['Jet Diameter (m):\t' num2str(D) '\n']);
    fprintf(fid,['Mach Number:\t' num2str(M) '\n']);
    fprintf(fid,['Ambient Temperature (K):\t' num2str(AT+273.15) '\n']);
    fprintf(fid,['Ambient Pressure (kPa):\t' num2str(AP) '\n']);
    fprintf(fid,['Relative Humidity (%%):\t' num2str(RH) '\n']);
    fprintf(fid,['Stagnation Temperature (K):\t' num2str(To) '\n']);
    fprintf(fid,['Exit Temperature (K):\t' num2str(Temp) '\n']);
    fprintf(fid,['Total Temperature Ratio:\t' num2str(To/(AT+273.15)) '\n']);
    fprintf(fid,['Exit Velocity (m/s):\t' num2str(M*sqrt(1.4*287.05*Temp)) '\n']);
    fprintf(fid,['Jet Velocity Ratio (Uj/a_inf):\t' num2str(M*sqrt(Temp/(AT+273.15))) '\n']);
    fprintf(fid,['Exit Density (kg/m^3):\t' num2str(1.01e5/(287.05*Temp)) '\n']);
    fprintf(fid,['Exit Viscosity (Pa s):\t' num2str((1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2))) '\n']);   %equation is Sutherlands formula - 291.15 is reference temperature
    fprintf(fid,['Jet Reynolds Number:\t' num2str(ReN,'%.0f') '\n']);
    
        %Calculates and prints OASPL for all channels    
    Bg = find(Std > 0.04,1,'first'); %Excludes spectral components below 0.04 Strouhal number
    Ed = find(Std > 4,1,'first');    %Excludes spectral components above 4 Strouhal number
    if isempty(Bg)
        Bg = 1;
    end
    if isempty(Ed)
        Ed = length(Std);
    end
    OASPL = 10*log10(trapz(fx(Bg:Ed),10.^(dBCal(Bg:Ed,:)/10),1)); clear dBCal;
    fprintf(fid,'OASPL (dB):'); fprintf(fid,'\t%.3f',OASPL); fprintf(fid,'\n');
        %Calculates and prints the average acoustic energy
    delta_Pol = abs(diff(ChPol)); IW = ([delta_Pol/2 0] + [0 delta_Pol/2])/(max(ChPol)-min(ChPol));
    AAE = 10*log10(sum(IW.*10.^(OASPL/10)));
    fprintf(fid,['Average Energy (dB):\t' num2str(AAE,'%.3f') '\n']);
        %Calculates and prints the normalized jet total power - assumes
        %azimuthally symmetric and fully expanded. The solid angle of the 
        %integration is divided out giving (W/sr). This power per solid
        %angle is then normalized by the jet mechanical power per solid
        %angle of a complete sphere (i.e. 4pi).
    NJP = 10*log10(32*NormD^2*pref^2*sqrt((AT+273.15)/Temp)*...
        trapz(ChPol*pi/180,10.^(OASPL/10).*sin(ChPol*pi/180))/...
        ((1.4*AP*1000)^2*M^3*(cos(ChPol(1)*pi/180)-cos(ChPol(end)*pi/180))));
    fprintf(fid,['Normalized Jet Power (dB):\t' num2str(NJP,'%.3f') '\n']);
    
        %Calculates and prints the peak SPL magnitude and Strouhal number
    fprintf(fid,'Peak of SPL Curve\n');    
    fprintf(fid,'  Magnitude (dB):');
    for nn = 1:Nch
        [c,I] = max(dBCal2(1:Ed,nn));
        I(I < 6) = 6;
        c = mean(dBCal2(I-5:I+5,nn));
        fprintf(fid,['\t' num2str(c,'%.3f')]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'  Frequency (Std):');
    for nn = 1:Nch
        [c,I] = max(dBCal2(1:Ed,nn));
        fprintf(fid,['\t' num2str(Std(I),'%.4f')]);
    end
    fprintf(fid,'\n');
        %Calculates and prints the high frequency slope of the spectra 
        %including uncertainty as seen on semilogx scale.
    fprintf(fid,'Slope of SPL Curve (Strouhal axis):');
    for nn = 1:Nch
        I = find(Std >= 0.5,1,'first');   %Ignores all points to left of Std = 0.7
        Pval = mmpolyfit(log10(Std(I:Ed)),dBCal2(I:Ed,nn),1,'Weight',1./Std(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
            %Calculates uncertainty of fit
        Serr = sum((Pval(1)*log10(Std(I:Ed))+Pval(2)-dBCal2(I:Ed,nn)).^2);
        merr = length(Std(I:Ed));
        Derr = merr*sum(log10(Std(I:Ed)).^2) - sum(log10(Std(I:Ed)))^2;
        PE = sqrt(Serr*merr/(merr-2)/Derr);
        
        fprintf(fid,['\t' num2str(Pval(1),'%10.3f') ' +/- ' num2str(PE,'%10.3f')]);
    end
    fprintf(fid,'\n');
    
        %If forcing is present, calculates/prints several relevant values
    if SHarmonics
        fprintf(fid,['Forcing Frequency (Hz):\t' num2str(FH) '\n']);
        fprintf(fid,['Forcing Strouhal Number:\t' num2str(FH*D/(M*sqrt(7/5*287.05*Temp))) '\n']);
            %Calculates and prints OASPL for all channels
        if DT~=0
            OASPL_DT = 10*log10(trapz(fx(Bg:Ed),10.^(dBCal2(Bg:Ed,:)/10),1)); clear dBCal2;
            fprintf(fid,'OASPL_Detoned (dB):'); fprintf(fid,'\t%.3f',OASPL_DT); fprintf(fid,'\n');
                %Calculates and prints the average acoustic energy
            AAE_DT = 10*log10(sum(IW.*10.^(OASPL_DT/10)));
            fprintf(fid,['Average Energy_Detoned (dB):\t' num2str(AAE_DT,'%.3f') '\n']);
                %Calculates and prints the normalized jet total power
            NJP_DT = 10*log10(32*NormD^2*pref^2*sqrt((AT+273.15)/Temp)*...
                trapz(ChPol*pi/180,10.^(OASPL_DT/10).*sin(ChPol*pi/180))/...
                ((1.4*AP*1000)^2*M^3*(cos(ChPol(1)*pi/180)-cos(ChPol(end)*pi/180))));
            fprintf(fid,['Normalized Jet Power_Detoned (dB):\t' num2str(NJP_DT,'%.3f') '\n']);     
        end
        
            %Finds SET file(s) for Baseline.  If none are found, the
            %following calculations are not performed.
        slist = struct2cell(dir(out_dir));
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
            OASPLB = zeros(1,Nch);  AAEB = zeros(1,length(slist));
            OASPLB_DT = OASPLB; AAEB_DT = AAEB; BT = AAEB; NJPB = AAEB; NJPB_DT = AAEB;
            
            if ~exist('SETpnts','var')
                    %Gets the SET file line numbers which contain the needed
                    %information.  
                q = strfind(slist{1},'_');
                if ~isempty(q)
                    fp{1} = slist{1}(1:q(1)-1);
                    if length(q) > 1
                        for nn = 2:length(q)
                            fp{nn} = slist{1}(q(nn-1)+1:q(nn)-1);
                        end
                    else
                        nn = 1;
                    end
                    fp{nn+1} = slist{1}(q(end)+1:end);

                    sOffset = 0;
                    for nn = 1:length(fp)
                        sOffset = sOffset+1;
                        if strmatch(fp{nn}(1),'M')   %Mach number is already in the SET file
                            sOffset = sOffset-1;
                        elseif strmatch(fp{nn}(1),'F')   %Forcing Frequency is already in the SET file
                            sOffset = sOffset-1;
                        elseif strmatch(fp{nn}(1),'T')   %Stagnation temperature is already in the SET file
                            sOffset = sOffset-1;
                        elseif strmatch(fp{nn}(1:2),'Ba')  %Ignores Baseline tag
                            sOffset = sOffset-1;
                        end
                    end
                end
                    
                done = false; SETpnts = zeros(1,7);
                fid2 = fopen([out_dir filesep slist{1}],'r');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(1) = SETpnts(1)+1;
                    done = ~isempty(strmatch('Stagnation Temperature (K):',L));
                end
                done = false; fseek(fid2,0,'bof');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(2) = SETpnts(2)+1;
                    done = ~isempty(strmatch('OASPL (dB):',L));
                end
                done = false; fseek(fid2,0,'bof');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(3) = SETpnts(3)+1;
                    done = ~isempty(strmatch('Average Energy (dB):',L));
                end
                done = false; fseek(fid2,0,'bof');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(4) = SETpnts(4)+1;
                    done = ~isempty(strmatch('Normalized Jet Power (dB):',L));
                end
                done = false; fseek(fid2,0,'bof');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(5) = SETpnts(5)+1;
                    done = ~isempty(strmatch('OASPL_Smoothed (dB):',L));
                end
                done = false; fseek(fid2,0,'bof');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(6) = SETpnts(6)+1;
                    done = ~isempty(strmatch('Average Energy_Smoothed (dB):',L));
                end
                done = false; fseek(fid2,0,'bof');
                while ~done
                    L = fgetl(fid2);
                    SETpnts(7) = SETpnts(7)+1;
                    done = ~isempty(strmatch('Normalized Jet Power_Smoothed (dB):',L));
                end
                fclose(fid2);
                SETpnts = SETpnts +sOffset;
            end
            
            for nn = 1:length(slist)    %reads baseline SET file(s) and extracts OASPL and total measured energy
                [NOASPLB,VB,SETpnts] = readSETfile([out_dir filesep slist{nn}],SETpnts);
                BT(nn) = VB{1}; %Baseline temperature (K)
                OASPLB(nn,:) = VB{2};
                AAEB(nn) = VB{3};
                NJPB(nn) = VB{4};
                OASPLB_DT(nn,:) = VB{5};
                AAEB_DT(nn) = VB{6};
                NJPB_DT(nn) = VB{7};
            end
            if length(slist) > 1
                [BT,IX] = sort(BT); %Ensures baseline data is in increasing order of temperature
                OASPLB = OASPLB(IX,:);
                OASPLB_DT = OASPLB_DT(IX,:);
                AAEB = AAEB(IX);
                AAEB_DT = AAEB_DT(IX);
                NJPB = NJPB(IX);
                NJPB_DT = NJPB_DT(IX);
                
                uBT = unique(BT);   %Averages data non-distinct temperatures
                if length(uBT)~=length(BT)
                    Tmp = cell(6,1);
                    for nn = 1:length(uBT)
                        q = find(BT==uBT(nn));
                        Tmp{1}(nn,:) = mean(OASPLB(q,:),1);
                        Tmp{2}(nn,:) = mean(OASPLB_DT(q,:),1);
                        Tmp{3}(nn) = mean(AAEB(q));
                        Tmp{4}(nn) = mean(AAEB_DT(q));
                        Tmp{5}(nn) = mean(NJPB(q));
                        Tmp{6}(nn) = mean(NJPB_DT(q));
                    end
                    BT = uBT;
                    OASPLB = Tmp{1};
                    OASPLB_DT = Tmp{2};
                    AAEB = Tmp{3};
                    AAEB_DT = Tmp{4};
                    NJPB = Tmp{5};
                    NJPB_DT = Tmp{6}; clear uBT Tmp;
                end
                
                for nn = 1:Nch      %Interpolate baseline results to most exactly match baseline characteristics at "To"
                    OASPLB(1,nn) = interp1(BT,OASPLB(:,nn),To,'linear','extrap');
                    OASPLB_DT(1,nn) = interp1(BT,OASPLB_DT(:,nn),To,'linear','extrap');
                end
                OASPLB = OASPLB(1,:);
                OASPLB_DT = OASPLB_DT(1,:);
                AAEB = interp1(BT,AAEB,To,'linear','extrap');
                AAEB_DT = interp1(BT,AAEB_DT,To,'linear','extrap');
                NJPB = interp1(BT,NJPB,To,'linear','extrap');
                NJPB_DT = interp1(BT,NJPB_DT,To,'linear','extrap');
            end
            
                %Calculates and prints delta OASPL
            fprintf(fid,'dOASPL (dB):'); fprintf(fid,'\t%.4f',OASPL - OASPLB); fprintf(fid,'\n');
                %Calculates and prints change in average acoustic energy
            fprintf(fid,['dAAE (dB):\t' num2str(AAE-AAEB,'%.4f') '\n']);
                %Calculates and prints change in normalized jet power
            fprintf(fid,['dNJP (dB):\t' num2str(NJP-NJPB,'%.4f') '\n']);
            
            if DT~=0
                    %Calculates and prints delta OASPL
                fprintf(fid,'dOASPL_Detoned (dB):'); fprintf(fid,'\t%.4f',OASPL_DT - OASPLB_DT); fprintf(fid,'\n');
                    %Calculates and prints change in average acoustic energy
                fprintf(fid,['dAAE_Detoned (dB):\t' num2str(AAE_DT-AAEB_DT,'%.4f') '\n']);
                    %Calculates and prints change in normalized jet power
                fprintf(fid,['dNJP_Detoned (dB):\t' num2str(NJP_DT-NJPB_DT,'%.4f') '\n']);
            end
        end       
    else
            %Calculates and prints OASPL for all channels    
        OASPL_DT = 10*log10(trapz(fx(Bg:Ed),10.^(dBCal2(Bg:Ed,:)/10),1)); clear dBCal;
        fprintf(fid,'OASPL_Smoothed (dB):'); fprintf(fid,'\t%.3f',OASPL_DT); fprintf(fid,'\n');
            %Calculates and prints the average acoustic energy
        AAE_DT = 10*log10(sum(IW.*10.^(OASPL_DT/10)));
        fprintf(fid,['Average Energy_Smoothed (dB):\t' num2str(AAE_DT,'%.3f') '\n']);
            %Calculates and prints the normalized jet total power
        NJP_DT = 10*log10(32*NormD^2*pref^2*sqrt((AT+273.15)/Temp)*...
            trapz(ChPol*pi/180,10.^(OASPL_DT/10).*sin(ChPol*pi/180))/...
            ((1.4*AP*1000)^2*M^3*(cos(ChPol(1)*pi/180)-cos(ChPol(end)*pi/180))));
        fprintf(fid,['Normalized Jet Power_Smoothed (dB):\t' num2str(NJP_DT,'%.3f') '\n']);
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECORD LOG INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('..Recording Log Information')
dep_dir = [out_dir filesep 'Dependencies'];
mkdir(dep_dir)

for n = 1:Nch   %Copies all calibration files to dependencies folder
    copyfile([CalPath filesep CalFiles{n}],[dep_dir filesep CalFiles{n}])
end

DEP = fdep(mfilename,'-q'); %Determines function dependency tree for program
for n = 1:length(DEP.fun)   %Copies all functions used to dependencies folder
    [q,nm,ext] = fileparts(DEP.fun{n});
    copyfile(DEP.fun{n},[dep_dir filesep nm ext]);
end
    %Writes log file containing list of files processed
fid = fopen([dep_dir filesep 'ProcessedFiles.txt'],'w');   
fprintf(fid,['Files processed on: ' datestr(tstamp,1) '\n']);
fprintf(fid,['Start Time: ' datestr(tstamp,13) ', End Time: ' datestr(now,13) '\n']);
for n = 1:length(flist)
    fprintf(fid,[flist{n} '\n']);
end
fclose(fid);

disp(' Done')   %Displays "Done" when all files have been processed

