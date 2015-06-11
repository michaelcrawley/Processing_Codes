function SoundProcess_v4()
%For the processing of acoustic data.  Program produces several files of
%processing results for an arbitrary number of channels of data.  Only one
%data directory may be processed per program execution.  Program requires
%initial (not automatable) user input to specify file locations and other
%parameters.  The program works on the assumption that there are N blocks
%of data in a file and that a RMS-FFT should be created from those N
%blocks.  It is assumed that any one directory being processed was
%collected at a constant ambient temperature.
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
%can be more than three characters).
% - If forcing harmonics are to be removed, file must have "_Fxxx" in the
% file name where xxx is a floating point value for the forcing frequency 
%in kHz (xxx can be more than three characters).
% - Baseline and forcing cases must be in the same directory. A baseline
% case should have "Baseline" in the file name.
% - File naming should follow this kind of convention
% "M0.9_m0_F010_T243.1_AUH.NOS" - The important thing is that each parameter be
% separated by an underscore. This example was for: Mach 0.9, mode zero
% forcing, 10 kHz forcing frequency, 243.1 Celcius jet stagnation
% temperature, and automatic duty cycle forcing.
% - A Corresponding .NOS calibration file for each channel must exist. 
%These files must all be in one directory, but it doesn't have to be the 
%same directory as the data.  These files must contain one column of data
%or N columns in ascending channel order where N is the number of channels.
% - All calibration files must have "_CHx" in there name where x is the
% number of the channel (x can be more than one character).
% - Calibration files should be noise data acquired when microphone is
% exposed to tone generator signal.  Program will assume that the peak
% value in the calibration spectrum is this reference signal.  The
% parameter ref is the amplitude (dB) of this signal.
% - Microphone distance measurements (x/D) must be known.  These will be
% manually inputted by the user during execution.
% - Ambient temperature (C) must be known.  This will be inputted by the
% user during execution.  It is assumed that a single value will suffice
% for all data being processed.
%
%    OUTPUT FILES
% - .fftNOS - File containing: Frequency axis, Strouhal number axis, and a
%column of SPL (dB) data for each channel acquired.
% - .HR.fftNOS - File containing: Frequency axis, Strouhal number axis, and
%a column of SPL (dB) data in which harmonics have been removed and profile 
%has been smoothed for each channel acquired.
% - .fig - MATLAB figure containing plots of SPL spectra for each channel
%acquired.
% - .png - standard picture file containing copy of image in .fig file.
% - .SET - Tab delimited file containing useful documented/calculated 
%parameters for the data.


    %Hard-coded parameters which may need to be changed
ref = 114;  %reference amplitude - dB
D = 0.0254; %jet diameter - m
M = 1.3;   %Mach Number
NormD = 80; %Normalization distance - x/D
NS = 1;     %0 = No tone subtraction
dBShift = 0;   %db Offset of spectra
ChPol = [90 80 70 60 50 45 40 35 30 25];   %Mic polar angles - degrees
ChD = [49 49.75 52.25 56.5 64 69.25 76.25 85.5 98 116]; %Mic distance - r/D

SD = 'F:\Research\Mach1.3\Acoustics\';   %Starting directory
	%Code assumes isentropically/ideally expanded flow

DPoints = 10;   %Trim first and last N points off spectra.
   
disp(' ')   %Asks user to verify proper parameters
disp(['Processing Running using Reference: ' num2str(ref) ' dB, D: ' num2str(D) ' m, NormD: ' num2str(NormD) ' x/D, and M: ' num2str(M)])
yn = input('  Is this Correct (y/n): ','s');
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
keep = zeros(size(q));
for n = 1:length(keep)
    if ~isempty(q{n})
        keep(n) = 1;
    end
end
flist = flist(1,logical(keep));
[keep,ok] = listdlg('PromptString','Select files to process:',...
                'SelectionMode','multiple',...
                'ListString',flist);    %Asks the user to select the subset of data files they wish to process
if ok==0
    error('PT_UserCancel','Program Terminated Due to User Selection of Cancel')
end
sampleRate = input('  Enter the Sample Rate for Data: ');   %Asks the user for the rate at which the data was acquired - If user enters nothing (Press Enter), program defaults to 200 kSamples/sec
if isempty(sampleRate)
    sampleRate = 200000;
    disp('    Sample Rate Set to Default: 200000')
end
dSize = input('  Enter the number of data points in a block: ');    %Asks the user for the number of data points in a block on which FFT will be calculated - If user enters nothing (Press Enter), program defaults to 8192 points
if isempty(dSize)
    dSize = 8192;
    disp('    Block Size Set to Default: 8192');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALIBRATION FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Asks user to select the calibration files to use and verifies proper
    %naming convention.
[FileName,PathName] = uigetfile('*.NOS','Select Calibration File(s) for Data',SD,'MultiSelect','on');
Nch = length(FileName); %Determines the number of channels of data by the number of calibration files
q = strfind(FileName, '_CH');
if isempty(q{1})
    error('Channel Number Prefix Must be All Caps...CH not Ch')
end
for n = 1:Nch   %Extracts channel number from from file name
    qq = FileName{n}(q{n}+3:q{n}+8);
    qqq = strfind(qq,'_');
    ChN(n) = str2num(qq(1:qqq(1)-1));
end
[B,IX] = sort(ChN);
FileName = FileName(IX);    %Reorders files in ascending channel number to ensure proper correlation to data columns

if ~(exist('ChPol','var') & exist('ChD','var'))
    for n = 1:Nch   %Asks user for distance to microphone for each channel
        ChPol(n) = input(['  Polar Angle (deg) for Channel ' num2str(n) ': ']);
        ChD(n) = input(['  Channel Distance (x/D) for Channel ' num2str(n) ': ']);
    end
elseif length(ChPol) ~= Nch
    disp('    Mic placement information mismatch!!')
    ChPol = []; ChD = [];
    for n = 1:Nch   %Asks user for distance to microphone for each channel
        ChPol(n) = input(['  Polar Angle (deg) for Channel ' num2str(n) ': ']);
        ChD(n) = input(['  Channel Distance (x/D) for Channel ' num2str(n) ': ']);
    end
end
for n = 1:Nch   %Calculates spectral reference value for each channel from the calibration files
    rRaw = dlmread([PathName FileName{n}],'\t');
    SrRaw = size(rRaw);
    if SrRaw(2) > 1 %If calibration file contains more than one column, program assumes nth column contains data for nth channel
        rRaw = rRaw(:,n);
    end
    nBlocks = length(rRaw)/dSize;   %Determines the number of blocks in calibration file data
    
        %Subtracts mean value of each block and applies a Hamming window
    rData = (rRaw(1:dSize)-mean(rRaw(1:dSize))).*hamming(dSize);  
    for m = 2:nBlocks
        rData(:,m) = (rRaw(dSize*(m-1)+1:dSize*m)-mean(rRaw(dSize*(m-1)+1:dSize*m))).*hamming(dSize);
    end
    rFFT = fft(rData,dSize,1);  %Calculates FFT
    rFFT = rFFT(1:ceil((dSize+1)/2),:); %Discards symmetric portion of spectrum
    mFFT = rFFT.*conj(rFFT)/dSize;  %Scales spectrum according to block size 
    if rem(dSize,2) %Conserves energy which would otherwise be lost by discardiing symmetric portion
        mFFT(2:end,:) = mFFT(2:end,:)*4;
    else
        mFFT(2:end-1,:) = mFFT(2:end-1,:)*4;
    end
    avgFFT = sqrt(sum(mFFT,2)/nBlocks); %Calculates RMS average of spectrum
    ChRef(n) = max(avgFFT); %Assumes peak of spectrum is reference value
end

fx = (0:ceil((dSize+1)/2)-1)*sampleRate/dSize;  %Calculates frequency axis
fx = fx(DPoints:end-DPoints); %First and last points will be thrown out as garbage 

color = 'brgkymc';  %Initializes array of colors for plotting

    %Inserts extra \ into directory path so it will print properly
Q = fliplr(strfind(PathName,'\'));
for nn = 1:length(Q)
    PathName = [PathName(1:Q(nn)) PathName(Q(nn):end)];
end
disp('  Calibration Data Processed...Beginning Data Processing')
AT = input('  Enter Ambient Temperature (C): ');    %Asks user for ambient temperature 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %All operations beyond this point are done without user interaction.
for n = 1:length(keep)
    Fnm = flist{keep(n)};
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
        FH = str2num(Fnm(FLoc+2:Fend-1))*1000;
    else
        SHarmonics=false;   
        FH = 0;
    end
    
    Temp = To/(1+1/5*M^2);  %Calculates isentropically expanded jet exit temperature - K
    Std = fx*D/(M*sqrt(7/5*287.05*Temp));    %Creates Strouhal number axis assuming Mach number is constant
    ReN = 1.01e5/(287.05*Temp) *M*sqrt(7/5*287.05*Temp)*...
        D /(1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2));  %Calculates jet exit Reynolds number using Sutherlands formula for viscosity and ideal gas law for density
    
    rRaw = dlmread([dir_name filesep flist{keep(n)}],'\t'); %reads data file
    S = size(rRaw);
    nBlocks = S(1)/dSize;   %Determines number of blocks in file - Note: All channels are processed simultaneously
    rData = rRaw(1:dSize,:)-repmat(mean(rRaw(1:dSize,:),1),[dSize,1]);  %Subtracts mean value from first block for every column
    rFFT = fft(rData,dSize,1);  %calculates FFT of first block for all channels 
    for m = 2:nBlocks
        rData = rRaw(dSize*(m-1)+1:dSize*m,:)-repmat(mean(rRaw(dSize*(m-1)+1:dSize*m,:),1),[dSize,1]);
        rFFT(:,:,m) = fft(rData,dSize,1);
    end
    rFFT = rFFT(1:ceil((dSize+1)/2),:,:);   %Discards symmetric portion
    mFFT = rFFT.*conj(rFFT)/dSize;  %Scales spectrum according to block size 
    if rem(dSize,2) %Conserves energy which would otherwise be lost by discardiing symmetric portion
        mFFT(2:end,:,:) = mFFT(2:end,:,:)*4;
    else
        mFFT(2:end-1,:,:) = mFFT(2:end-1,:,:)*4;
    end
    avgFFT = sqrt(sum(mFFT,3)/nBlocks); %Calculates RMS average spectrum for all channels
    dB = ref+ 20*log10(avgFFT./repmat(ChRef,[ceil((dSize+1)/2),1]));    %Converts into dB using calibration and reference values
    dBCal = dB + 20*log10(repmat(ChD,[ceil((dSize+1)/2),1])/NormD) +dBShift;     %Normalizes data to a standard distance

    dBCal = dBCal(DPoints:end-DPoints,:);   %removes junk points at beginning and end
      
    dBCal2 = dBCal;
    if NS~=0
        for m = 1:Nch   %Smoothes and removes forcing harmonics from spectrum. If FH==0 (forcing frequency) program doesn't look for harmonics
            dBCal2(:,m) = ToneSubtract_v3(fx',dBCal(:,m),FH,10);
        end
    else
        for m = 1:Nch   %Smoothes and removes forcing harmonics from spectrum. If FH==0 (forcing frequency) program doesn't look for harmonics
            dBCal2(:,m) = ToneSubtract_v3(fx',dBCal(:,m),0,10);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OUTPUT FILE GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%  PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure  %Plots all spectra on a single figure
    if Nch > 7
        semilogx(Std,dBCal(:,1),'Color',rand(3,1))
        LBL = 'Channel 1';  %Variable contains labels which will become column headers on output files
        hold on
        for m = 2:Nch
            semilogx(Std,dBCal(:,m),'Color',rand(3,1))
            LBL = [LBL '\tChannel ' num2str(m)];
        end
    else
        semilogx(Std,dBCal(:,1),color(1))
        LBL = 'Channel 1';
        hold on
        for m = 2:Nch
            semilogx(Std,dBCal(:,m),color(m))
            LBL = [LBL '\tChannel ' num2str(m)];
        end
    end
    for m = 1:Nch   %Plots smoothed spectra on top of original spectra
        semilogx(Std,dBCal2(:,m),color(mod(m,7)+1))
        LEGN{m} = num2str(m);
    end
    legend(LEGN)
    xlabel(['Strouhal Number for D = ' num2str(D) ' m Jet at M = ' num2str(M)])
    ylabel('SPL (dB)')
    title({['Average Spectrum for: ' flist{keep(n)}],['Data Normalized for ' num2str(NormD) ' x/D'],['Reynolds Number: ' num2str(ReN)]})
    grid on
    saveas(gcf,[dir_name filesep flist{keep(n)}(1:end-4) '_Std.png'])   %Saves figure as .png
    saveas(gcf,[dir_name filesep flist{keep(n)}(1:end-4) '_Std.fig'])   %Saves figure as .fig
    close    %Closes figure
    
    %%%%%  SPECTRAL DATA FILES  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([dir_name filesep flist{keep(n)}(1:end-3) 'fftNOS'],'w');   %Opens .fftNOS file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number\t' LBL '\n']);
    fclose(fid);
    dlmwrite([dir_name filesep flist{keep(n)}(1:end-3) 'fftNOS'],[fx' Std' dBCal], 'delimiter', '\t', '-append');   %writes spectral data to file
    clear dbCal %All SET file calculations use the smoothed spectra
    
    fid = fopen([dir_name filesep flist{keep(n)}(1:end-3) 'HR.fftNOS'],'w');   %Opens .HR.fftNOS file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number\t' LBL '\n']);
    fclose(fid);
    dlmwrite([dir_name filesep flist{keep(n)}(1:end-3) 'HR.fftNOS'],[fx' Std' dBCal2], 'delimiter', '\t', '-append');   %writes smoothed spectral data to file

    %%%%%  SET FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([dir_name filesep flist{keep(n)}(1:end-3) 'SET'],'w');
    fprintf(fid,'PROCESS VALUES FOR DATA\n');
    fprintf(fid,['Sampling Rate (Hz):\t' num2str(sampleRate) '\n']);
    fprintf(fid,['Block Size:\t' num2str(dSize) '\n']);
    fprintf(fid,['Reference Amplitude (dB):\t' num2str(ref) '\n']);
    fprintf(fid,['Jet Diameter (m):\t' num2str(D) '\n']);
    fprintf(fid,['Mach Number:\t' num2str(M) '\n']);
    fprintf(fid,['Normalization Distance (x/D):\t' num2str(NormD) '\n']);    
    fprintf(fid,['Number of Channels:\t' num2str(Nch) '\n']);
        %Records location and name of each calibration file used.
    fprintf(fid,['Location of Calibration Files:\t' PathName '\n']);
    fprintf(fid,'Calibration Files:\t');
    for nn = 1:length(FileName)
        fprintf(fid,[FileName{nn} '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Channel Distance (x/D):\t');
    for nn = 1:length(ChD)
        fprintf(fid,[num2str(ChD(nn)) '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Channel Polar Angle (deg):\t');
    for nn = 1:length(ChPol)
        fprintf(fid,[num2str(ChPol(nn)) '\t']);
    end
    fprintf(fid,'\n');
    
    fprintf(fid,['Stagnation Temperature (K):\t' num2str(To) '\n']);
    fprintf(fid,['Ambient Temperature (K):\t' num2str(AT+273.15) '\n']);
    fprintf(fid,['Jet Temperature (K):\t' num2str(Temp) '\n']);
    fprintf(fid,['Jet Temperature Ratio:\t' num2str(To/(AT+273.15)) '\n']);
    fprintf(fid,['Jet Velocity (m/s):\t' num2str(M*sqrt(7/5*287.05*Temp)) '\n']);
    fprintf(fid,['Jet Velocity Ratio:\t' num2str(M*sqrt(Temp/(AT+273.15))) '\n']);
    fprintf(fid,['Jet Density (kg/m^3):\t' num2str(1.01e5/(287.05*Temp)) '\n']);
    fprintf(fid,['Jet Viscosity (Pa s):\t' num2str((1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2))) '\n']);   %equation is Sutherlands formula - 291.15 is reference temperature
    fprintf(fid,['Reynolds Number:\t' num2str(ReN) '\n']);
        %Calculates and prints OASPL for all channels
    fprintf(fid,['OASPL:\t']);      
    Bg = min(find(Std > 0.04)); %Excludes spectral components below 0.04 Strouhal number
    Ed = min(find(Std > 4));    %Excludes spectral components above 4 Strouhal number
    for nn = 1:length(ChD)
        fprintf(fid,[num2str(10*log10(sum(10.^(dBCal2(Bg:Ed,nn)/10)))) '\t']);
    end
    fprintf(fid,'\n');
        %Calculates and prints the average acoustic energy
    fprintf(fid,['Average Energy (dB):\t']);       
    delta_Pol = abs(diff(ChPol)); IW = ([delta_Pol/2 0] + [0 delta_Pol/2])/(max(ChPol)-min(ChPol)); 
    fprintf(fid,[num2str(10*log10(sum(IW.*sum(10.^(dBCal2(Bg:Ed,:)/10))))) '\n']);
        %Calculates and prints the peak SPL magnitude and Strouhal number
    fprintf(fid,'Peak of SPL Curve\n');
    fprintf(fid,['  Magnitude (dB):\t']);
    for nn = 1:length(ChD)
        [c,I] = max(dBCal2(1:Ed,nn));
        I(I < 6) = 6;
        c = mean(dBCal2(I-5:I+5,nn));
        fprintf(fid,[num2str(c) '\t']);
    end
    fprintf(fid,'\n');
    fprintf(fid,['  Frequency (Std):\t']);
    for nn = 1:length(ChD)
        [c,I] = max(dBCal2(1:Ed,nn));
        fprintf(fid,[num2str(Std(I)) '\t']);
    end
    fprintf(fid,'\n');
        %Calculates and prints the high frequency slope of the spectra 
        %including uncertainty as seen on semilogx scale.
    fprintf(fid,['Slope of SPL Curve (Strouhal axis):\t']);
    for nn = 1:length(ChD)
        [c,I] = max(dBCal2(1:Ed,nn));   %Ignores all points to left of peak
        Pval = mmpolyfit(log(Std(I*2:Ed)),dBCal2(I*2:Ed,nn),1,'Weight',1./Std(I*2:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
            %Calculates uncertainty of fit
        Serr = sum((Pval(1)*log(Std(I*2:Ed)')+Pval(2)-dBCal2(I*2:Ed,nn)).^2);
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
            %Finds SET file(s) for Baseline.  If none are found, the
            %following calculations are not performed.
        slist = struct2cell(dir(dir_name));
        qs = strfind(slist(1,:),'Baseline');
        qs2 = strfind(slist(1,:),'.SET');
        kp = zeros(size(qs));
        for n = 1:length(qs)
            if (~isempty(qs{n}))&&(~isempty(qs2{n}))
                kp(n) = 1;
            end
        end
        if sum(kp) > 0  %If baseline SET file(s) are found
            slist = slist(1,logical(kp));
            OASPLB = [];
            AAE = [];
            for nn = 1:length(slist)    %reads baseline SET file(s) and extracts OASPL and total measured energy
                [NOASPLB,VB,kp] = readSETfile([dir_name filesep slist{nn}],[21 22]);
                OASPLB = [OASPLB; VB{1}];
                AAE = [AAE; VB{2}];
            end
            OASPLB = mean(OASPLB,1);    %averages OASPL from multiple baselines
            AAE = mean(AAE);  %averages average acoustic energy from multiple baselines
            
                %Calculates and prints delta OASPL
            fprintf(fid,['dOASPL:\t']);      
            for nn = 1:length(ChD)
                fprintf(fid,[num2str(10*log10(sum(10.^(dBCal2(Bg:Ed,nn)/10))) - OASPLB(nn)) '\t']);
            end
            fprintf(fid,'\n');
                %Calculates and prints change in average acoustic energy
            fprintf(fid,['dAAE:\t']);      
            fprintf(fid,[num2str(10*log10(sum(IW.*sum(10.^(dBCal2(Bg:Ed,:)/10))))-AAE) '\n']);
        end        
    end
    fclose(fid);
end
disp(' Done')   %Displays "Done" when all files have been processed

