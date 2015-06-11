function SoundProcess_v14(PARAMS)
%For the processing of acoustic data.  Program produces several files of
%processing results for an arbitrary number of channels of data.  Only one
%data directory may be processed per program execution.  Program requires
%complete information in the PARAMS .mat file or initial user input to 
%specify file locations and other parameters. With a fully specified PARAMS 
%file, the program can operate with no user intervention. The program works 
%on the assumption that there are N blocks of data in a file and that a 
%RMS-FFT should be created from those N blocks.  It is assumed that any one 
%directory being processed was collected at a constant ambient temperature. 
%Code assumes isentropically/ideally expanded flow.
%
%INPUTS
%  PARAMS - The path to the .mat file containing the parameters. This file 
%   must always be used, but some of the variables in it are optional (see 
%   below).
%
%OUTPUT FILES
% - .fft* - File containing: Frequency axis, Strouhal number axis, and a
% column of SPL (dB) data for each channel acquired.
% - .S.fft* - File containing: Frequency axis, Strouhal number axis, and
% a column of SPL (dB) data in which harmonics have been removed and profile 
% has been smoothed for each channel acquired.
% - .fig - MATLAB figure containing plots of SPL spectra for each channel
% acquired.
% - .png - standard picture file containing copy of image in .fig file.
% - .SET - Tab delimited file containing useful documented/calculated 
% parameters for the data.
%
%
%%%%% PARAMS File Requirements/Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FLOW PROPERTY INPUTS
% D: jet diameter - m
% M: Hydrodynamic Mach Number
% AP: Atmospheric Pressure - kPa
% AT: Ambient Temperature - C
% RH: Relative Humidity - %
%
%DAQ PARAMETERS
% sampleRate: data sampling rate - Hz
% BS: Samples per block
% ChD: Microphone radial distance from nozzle exit - m
% ChPol: Microphone polar angles relative to jet axis - degrees
% MO: Mic orientation positive angle relative to incident sound - degrees
%       NOTE- MO==0 is normal incidence. Enter zero for no microphone 
%       orientation correction.
%
%FILE SELECTION
% SD: Starting directory for file selection. Unnecessary if data and cal
%   files are specified.
% EXT: File extension for data and cal files
%%%OPTIONAL VARIABLES%%%
% src_dir: The path to the directory containing the data files. This
%   variable is optional. If it doesn't exist in the 'PARAMS' file, the
%   program will ask the user to provide the path.
% flist: A cell array containing the list of files to be processed. This
%   variable is optional. If it doesn't exist in the 'PARAMS' file, the
%   program will ask the user to provide the list of files.
% CalPath: The path to the directory containing the calibration files. This
%   variable is optional. If it doesn't exist in the 'PARAMS' file, the
%   program will ask the user to provide the path.
% CalFiles: A cell array containing the list of calibration files. This
%   variable is optional. If it doesn't exist in the 'PARAMS' file, the
%   program will ask the user to provide the list of files.
%%%%%%%%%%%%%%%%%%%%%%%%
% ref: Reference amplitude of microphone calibration - dB
%
%SPECTRUM CALCULATION/CORRECTION
% NDs: Normalization distance correction method 
%   Accepted Inputs:'Radial'- simple radial correction
%                   'NS'- noise source location correction
%                   'None'- perform no propagation distance correction
% NormD: Normalization distance - x/D
% WNDO: Data windowing method. Use "rectwin" for no window.
% DPoints: Trim first and last N points off spectra
% DT: Logical switch for performing detoning. 0==No tone subtraction.
% atmAbs: Atmospheric absorption correction 
%   Accepted Inputs:'Standard'- correct to standard day
%                   'Lossless'- removes all absorption
%                   'None'- perform no absorption correction
% dBShift: A constant spectrum offset - dB
%       NOTE- This arbitrary offset can be used to correct spectra for
%       uncalibrated amplification levels. E.g. If data acquired at 100
%       mV/Pa and calibration data at 316 mV/Pa, set 'dBShift' to 10. It
%       can be a single value or a vector with a unique value for each
%       channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GENERAL PROGRAM REQUIREMENTS
% - All files must have an extension given in "EXT" in order to be processed.
% - All files must be tab delimited and contain only data.
% - All data for a given set of parameters must be in a single file (i.e.
% no spreading channels among multiple files).
% - Data must be arranged so that each column represents a channel.
% - Data must have an integer number of blocks (i.e. if there are 1000
% samples per block and 10 blocks, there better be 10000 rows in the file -
% no more, no less).  One consequence of this is that all channels must
% have the same sampling rate, samples per block, and number of blocks.
% - All data files must have "_Txxx" in the file name where xxx is a
% floating point value for the jet stagnation temperature in Celcius (xxx 
% can be more than three characters).
% - If forcing harmonics are to be removed, file must have "_Fxxx" in the
% file name where xxx is a floating point value for the forcing frequency 
% in kHz (xxx can be more than three characters).
% - Baseline and forcing cases must be in the same directory. A baseline
% case should have the following form: "Mxx_Baseline_Txxx.*".
% - File naming should follow this kind of convention
% "M0.9_m0_F010_T243.1_DAUH.*" - The important thing is that each parameter be
% separated by an underscore. This example was for: Mach 0.9, mode zero
% forcing, 10 kHz forcing frequency, 243.1 Celcius jet stagnation
% temperature, and automatic duty cycle forcing.
% - A Corresponding "EXT" calibration file for each channel must exist. 
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VERSION/UPDATE INFORMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-v14
%  - Changed calibration file reading process to handle binary files with
%  one or more channels.
%-v13
%  - Changed dbShift to accept a vector for correcting different
%  uncalibrated gains on different channels. This required an update to the
%  subroutine "NormDistance".
%  - Changed plotting color scheme to equally divide the jet colormap among
%  the number of channels.
%  - Improved subroutine "IDBaselines" to accept Strouhal number of zero
%  and to strip forcing parameters while leaving suffixes untouched.
%-v12
%  - Changed microphone distance variable "ChD" so that is inputted in
%  meters so that multiple nozzle diameters can be easily accomodated
%  without the need to change all the microphone distance values.
%-v11b
%  - Switched user confirmation of parameters to a dialog box. New
%  implementation uses a modified version of "questdlg" called
%  "questdlg_timeout" that has a timeout on the dialog box. If there is no
%  user input after 60 seconds, the program assumes that the function call
%  to "SoundProcess" is for a full-auto environment and proceeds with 
%  processing. This allows "SoundProcess" to be called in a for loop using
%  a set of fully specified PARAMS files. Since no user interaction is 
%  required, multiple sets of data can be processed in a single user
%  interaction.
%-v11
%  - Added ability to read binary data files. Assumes binary data files are
%  float32 and organized as one block for all channels then the next block
%  (e.g. 8192 samples for ch1, 8192 samples for ch2,...8192 samples for 
%  chN, next 8192 samples for ch1...).
%  - Added ability to read binary calibration files. Assumes that
%  calibration files contain only one channel of data and no timestamp.
%  - Program can function with binary data files and delimited calibration
%  files or vice-versa.
%  - Added a hard-coded parameter to strip out timestamp information if it
%  present.
%  - Added the ability to parse forcing Strouhal numbers "_S" in the file
%  name.
%-v10
%  - Added additional option of no correction for atmospheric absorption
%  - Reduced the computational overhead of calibration data by computing
%  prms^2 directly as opposed to computing OASPL from spectrum.
%  - Added additional option of no correction for propagation distance.
%  Note that this option will cause the reported normalized distance and 
%  jet power to return NaN since these values have no meaning in this
%  context.
%  - Added additional option of no correction for microphone orientation.
%  - Added parameter "EXT" for processing data with an arbitrary extension.
%  - Changed settings configuration to work from a "PARAMS.mat" file. File
%  can optionally include directory and file lists for calibration files
%  and/or data files. If not included, the program will prompt the user to
%  specify those files.
%-v9
%  - Saves all results into a directory whose name is a date-time stamp 
%  created at run-time.
%  - Produces a function dependency report and creates an archive 
%  containing: all m-files used to process the data, the calibration files, 
%  and a list of which files were processed. 
%  - Baseline data is now linearly interpolated versus temperature to more 
%  accurately match the characteristics at any given forced case.
%  - Smoothed spectra are now computed as a 30th order polynomial fit which
%  works better than moving average detoning for narrowband peaks (see
%  ToneSubtract_v5.m).
%  - Computes the Normalized Jet Power (NJP) and writes that additional
%  information to SET files.
%  - Computes correction for atmospheric with switch for choosing lossless
%  propagation or standard day propagation.
%  - Computes normalization distance propagation with switch for
%  propagation correction which accounts for a noise source distribution.
%  See function "NormDistance.m".
%  - Dynamically scans baseline SET files for locations of needed
%  information instead of assuming known locations.
%  - Properly ignores additional filename parts (e.g. Mxx_EXAMPLE_xxxx.*)
%  so that additional file identifiers may be used.
%-v8 
%  - will accept directories containing baselines whose file name uses a
%  forcing frequency of value zero. The original files with form
%  (Mxx_xxxx_F0.0_xxxx_Txxxx_xxxx.*) will be renamed as
%  (Mxx_Baseline_Txxxx.*) before the file list is presented to the user.
%-v7 
%  - modifies OASPL calculations to use non-smoothed spectra as much as 
%  possible and to compute dOASPL for both detoned and non-detoned forced
%  spectra if user requested detoning.
%-v6 
%  - adds the correction for atmospheric absorption. It updates the first 
%  user call to verify hard-coded parameters. It corrects the method for
%  calculating the FFT. Previous versions were calculating SPL levels which
%  were dependent on sampling properties.
%-v5 
%  - adds the correction for microphone orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChPol = []; clear ChPol;    %This line is necessary to fool MATLAB's error checking. Program will still run fine if commented, but the editor will tell you that there is a problem.
load(PARAMS); clear USAGE;

    %Examples of parameters to be submitted in "PARAMS" file which may need to be changed
% M = 1.3;   %Mach Number
% AT = 28;    %Ambient Temperature - C
% AP = 100.5; %Atmospheric Pressure - kPa
% RH = 50;    %Relative Humidity - %
% D = 0.0254; %jet diameter - m
% ref = 114;  %reference amplitude - dB
% sampleRate = 200000;    %data sampling rate - Hz
% BS = 8192;   %Samples per block
% NDs = 'Radial';     %Normalization distance correction method - inputs: 'Radial'-simple radial correction; 'NS'-noise source location correction; 'None'-perform no propagation distance correction
% NormD = 80; %Normalization distance - x/D
% ChPol = [90 80 70 60 50 45 40 35 30 25];   %Mic polar angles - degrees
% MO = [0 10 20 30 40 45 50 55 60 0];     %Mic orientation positive angle relative to incident sound - degrees - Note: MO==0 is normal incidence. Enter zero for no microphone orientation correction
% ChD = [49 49.75 52.25 56.5 64 69.25 76.25 85.5 98 116]; %Mic distance - m
% DT = 1;     %0 = No tone subtraction (DeToning)
% dBShift = 0;   %db Offset of spectra - this arbitrary offset can be used to correct spectra for uncalibrated amplification levels
% DPoints = 6;   %Trim first and last N points off spectra.
% WNDO = 'hamming';    %Data windowing method - use "rectwin" for no window
% atmAbs = 'Standard';    %Atmospheric absorption correction - inputs: 'Standard'-correct to standard day; 'Lossless'-removes absorption; 'None'-perform no absorption correction
% SD = 'F:\Research\';   %Starting directory
% EXT = 'NOS';    %File extension for data files

TS = false;		%false == no timestamp present in data files
TS_cal = false;	%false == no timestamp present in calibration files
delimI = '\t';	%Delimiter to check for in file reading 
delimO = '\t';	%Delimiter for file writing - As of 20101006 other functions assume output files are tab delimited...change at own risk.
fLimits = [0.01 4]; %The limits of the plotting on the figure windows. Use [xx xx] to limit only abscissa. Use [xx xx xx xx] to limit all dimensions.
iLimits = [0.04 4]; %The limits of integration (in Strouhal number) for determining OASPL.
pref = 20e-6;   %Pascals - reference pressure - used in calculating normalized jet power

if strcmpi('None',NDs)  %Forces Normalization distance to invalid if no propagation correction is wanted
    NormD = NaN;
end
if exist('src_dir','var') && exist('flist','var') && exist('CalPath','var') && exist('CalFiles','var')
    SD = 'N.A.';
end

    %Asks user to verify proper parameters
disp('Current Processing Values: ');
disp(['                 Ambient Temp (C): ' num2str(AT)]);
disp(['           Ambient Pressure (kPa): ' num2str(AP)]);
disp(['            Relative Humidity (%): ' num2str(RH)]);
disp(['                 Jet Diameter (m): ' num2str(D)]);
disp(['                      Mach Number: ' num2str(M)]);
disp(['      Reference Signal Level (dB): ' num2str(ref)]);
disp(['               Sampling Rate (Hz): ' num2str(sampleRate)]);
disp(['                Samples Per Block: ' num2str(BS)]);
disp(['Mic Normalization Distance Method: ' NDs]);
disp([' Mic Normalization Distance (r/D): ' num2str(NormD)]);
disp(['       Mic Polar Angles (degrees): ' num2str(ChPol,'%.0f ')]);
disp(['  Mic Orientation Angle (degrees): ' num2str(MO,'%.0f ')]);
disp(['                 Mic Distance (m): ' num2str(ChD,'%.2f ')]);
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
if exist('src_dir','var') 
    disp(['            Data Source Directory: ' src_dir]); %#ok<NODEF>
end
if exist('flist','var')
    disp('                   Data File List: PRESENT');
else
    disp('                   Data File List: NOT PRESENT');
end
if exist('CalPath','var')
    disp(['       Calibration File Directory: ' CalPath]); %#ok<NODEF>
end
if exist('CalFiles','var')
    disp('            Calibration File List: PRESENT');
else
    disp('            Calibration File List: NOT PRESENT');
end

yn = questdlg_timeout('Is the information correct? - Program assumes yes if no user selection after 60 seconds.','Parameter Confirmation',60);
if ~strcmpi(yn,'Yes') && ~strcmpi(yn,'')
    error('Fix Values in code and run again')
end
if strcmpi(yn,'')
	disp('  NO user confirmation of parameters...Assuming full-auto operation')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FILE SELECTION/BASIC CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Asks user for a directory containing data - can only process one
    %directory per program execution.
if ~exist('src_dir','var')
    src_dir = uigetdir(SD,'Specify Directory Containing Data'); 
end
USelect = false;    %switch variable that allows user to select a file list subset if the file list was not submitted in the "PARAMS" file
if ~exist('flist','var')
    USelect = true;
    flist = struct2cell(dir(src_dir));
    tmp = strfind(flist(1,:),['.' EXT]); %extracts list of data files - ignores all other files
    keep = ~cellfun('isempty',tmp); clear tmp;
    flist = flist(1,logical(keep));
end
flist = IDBaselines(src_dir,flist);    %rewrites file names of Baselines
if USelect
    [keep,ok] = listdlg('PromptString','Select files to process:',...
                    'SelectionMode','multiple',...
                    'ListString',flist);    %Asks the user to select the subset of data files they wish to process
    if ok==0
        error('Program Terminated Due to User Selection of Cancel')
    end
    flist = flist(keep); clear keep ok;
end
clear USelect;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALIBRATION FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('CalPath','var') || ~exist('CalFiles','var')
    %Asks user to select the calibration files to use and verifies proper
    %naming convention.
    [CalFiles,CalPath] = uigetfile(['*.' EXT],'Select Calibration File(s) for Data',SD,'MultiSelect','on');
end
%%%%% All operations beyond this point are done without user interaction. %%%%%

Nch = length(CalFiles); %Determines the number of channels of data by the number of calibration files
ChN = zeros(1,Nch);
for n = 1:Nch   %Extracts channel number from from file name
    tmp = regexpi(CalFiles{n},'_CH');
	if isempty(tmp)
        error('Channel Number Prefix Must be CH, Ch, or ch')
	end
    tmp2 = CalFiles{n}(tmp+3:tmp+8);
    tmp3 = strfind(tmp2,'_');
    ChN(n) = str2double(tmp2(1:tmp3(1)-1));
end
if length(unique(ChN))~=length(ChN)
    error('Calibration channel numbers are not unique')
end
[B,IX] = sort(ChN); clear n ChN B tmp tmp2 tmp3;
CalFiles = CalFiles(IX); clear IX; %Reorders files in ascending channel number to ensure proper correlation to data columns

if length(ChPol)~=Nch
    error('    Microphone information mismatch!! Number of calibration files does not match microphone location information.')
end

ChRef = zeros(1,Nch);
for n = 1:Nch   %Calculates spectral reference value for each channel from the calibration files
	if ~exist('isbinary','var')
		isbinary = true;
		fid = fopen([CalPath filesep CalFiles{n}],'r');
		tmp = fread(fid,500,'*char'); fclose(fid);
		if strcmp('\t',delimI)
			if length(strfind(tmp',char(9))) > 5
				isbinary = false;
			end
		else
			if length(strfind(tmp',delimI)) > 5
				isbinary = false;
			end
		end
		if isbinary
			if length(strfind(tmp',char(10))) > 5	%checks for line feed characters if delimiter wasn't found
				isbinary = false;
			end
		end
	end
	if isbinary
		fid = fopen([CalPath filesep CalFiles{n}],'r');
		tmp = fread(fid,'float32'); fclose(fid);
		
		Ltmp = length(tmp);
		if round(Ltmp/BS/Nch)==Ltmp/BS/Nch	%If there could be more than one channel, unpack as if there are "Nch" channels
			if TS_cal	%If timestamp is present in files
				tmp2 = reshape(tmp,BS,Nch+1,[]);	%parses data into channels along second dimension and blocks along third dimension
				tmp2 = permute(tmp2,[1 3 2]);	%reorders data into channels along third and blocks along second dimension
				tmp2 = reshape(tmp2,[],Nch+1);	%reshapes data into 2-D matrix with all blocks for a channel in one column
				tmp2 = tmp2(:,2:end);	%removes timestamp
			else
				tmp2 = reshape(tmp,BS,Nch,[]);
				tmp2 = permute(tmp2,[1 3 2]);
				tmp2 = reshape(tmp2,[],Nch);
			end
			mtmp = mean(tmp2.^2);	%If it was actually only one channel, all imagined channels will have same mean square. In which case the sum below will equal the number of channles. 
			if sum(round(mtmp/mtmp(n))) ~= Nch
				tmp = tmp2(:,n);
				disp(['     Cal File ' num2str(n) ' has ' num2str(Nch) ' channels in it']);
			end
			clear tmp2 mtmp;
		else
			disp(['     Cal File ' num2str(n) ' has 1 channel in it']);
		end
	else
		tmp = dlmread([CalPath filesep CalFiles{n}],delimI);
		if size(tmp,2)~=1   %If there is only one channel of data, assume it is for the nth channel    
			%Else, assume that the nth column contains the calibration for the nth channel
			tmp = tmp(:,n);
		end
	end
	S = length(tmp);
	if round(S/BS)~=S/BS
		error(['Data in ' CalFiles{n} ' is not an integer number of blocks'])
	end
    tmp = reshape(tmp,BS,[]);    %Cut the data into blocks
    tmp = std(tmp,0,1);     %Remove mean and compute RMS of blocks
    ChRef(n) = mean(tmp.^2);    %Compute mean square of prms of blocks (i.e. prms^2)
end
clear n S tmp isbinary;
ChRef = spdiags(1./ChRef(:),0,Nch,Nch); %prepares the calibration values for scaling the data.
fx = (0:ceil((BS+1)/2)-1)'*(sampleRate/BS);   %Frequency axis
fx = fx(DPoints:end-DPoints); %First and last points will be thrown out as garbage

    %Inserts extra \ into directory path so it will print properly
tmp = fliplr(strfind(CalPath,'\')); CalPathP = CalPath;
for n = 1:length(tmp)
    CalPathP = [CalPathP(1:tmp(n)) CalPathP(tmp(n):end)];
end
clear n tmp;
disp('  Calibration Data Processed...Beginning Data Processing')

tstamp = now;
out_dir = [src_dir filesep datestr(tstamp,30)];
mkdir(out_dir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:length(flist)
    Fnm = flist{n};
    disp(['    Processing File: ' Fnm])  %Displays file currently being processed
        %Locates temperature information in file name
    TLoc = strfind(Fnm,'_T');   
    Tend = min(strfind(Fnm(TLoc+1:end),'_'))+TLoc;
    if isempty(Tend)
        Tend = length(Fnm)-length(['.' EXT])+1;
    end
    Ttxt = Fnm(TLoc+2:Tend-1); clear TLoc Tend;
    Ttxt(strfind(Ttxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
    To = str2num(Ttxt)+273.15; clear Ttxt; %#ok<ST2NM> %Converts to Kelvin
	Temp = To/(1+1/5*M^2);  %Calculates isentropically expanded jet exit temperature - K
	
		%Locates forcing frequency information in file name if it exists
	FLoc = strfind(Fnm,'_F');  
	SLoc = strfind(Fnm,'_S');
    if ~isempty(FLoc)	%Looks for Forcing frequency (kHz) in filename
        SHarmonics=true;    %This boolean will be used to determine if forcing frequency information needs to be written in output file
        Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
        FH = str2num(Fnm(FLoc+2:Fend-1))*1000; %#ok<ST2NM>
		StDF = FH*D/(M*sqrt(7/5*287.05*Temp));
	elseif ~isempty(SLoc)	%Looks for Forcing strouhal number in filename
        SHarmonics=true;    %This boolean will be used to determine if forcing frequency information needs to be written in output file
        Fend = min(strfind(Fnm(SLoc+1:end),'_'))+SLoc;
        StDF = str2num(Fnm(SLoc+2:Fend-1)); %#ok<ST2NM>
		FH = StDF*M*sqrt(7/5*287.05*Temp)/D;
	else
        SHarmonics=false;   
        FH = 0;
    end
    clear FLoc SLoc Fend;
    
    Std = fx*D/(M*sqrt(7/5*287.05*Temp));    %Creates Strouhal number axis assuming Mach number is constant
    ReN = 1.01e5/(287.05*Temp) *M*sqrt(7/5*287.05*Temp)*...
        D /(1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2));  %Calculates jet exit Reynolds number using Sutherlands formula for viscosity and ideal gas law for density
    
	if ~exist('isbinary','var')	%only tests for binary file content based on first file read
		isbinary = true;
		fid = fopen([src_dir filesep Fnm],'r');
		rRaw = fread(fid,500,'*char'); fclose(fid);	%Reads first 200 bytes assuming they are ascii
		if strcmp('\t',delimI)
			if length(strfind(rRaw',char(9))) > 5	%If delimiter is present, file is not binary
				isbinary = false;
			end
		else
			if length(strfind(rRaw',delimI)) > 5
				isbinary = false;
			end
		end
		if isbinary
			if length(strfind(rRaw',char(10))) > 5
				isbinary = false;
			end
		end
	end
	if isbinary
		fid = fopen([src_dir filesep Fnm],'r');
		rRaw = fread(fid,'float32'); fclose(fid);
		if TS	%If timestamp is present in files
			rRaw = reshape(rRaw,BS,Nch+1,[]);	%parses data into channels along second dimension and blocks along third dimension
			rRaw = permute(rRaw,[1 3 2]);	%reorders data into channels along third and blocks along second dimension
			rRaw = reshape(rRaw,[],Nch+1);	%reshapes data into 2-D matrix with all blocks for a channel in one column
			rRaw = rRaw(:,2:end);	%removes timestamp
		else
			rRaw = reshape(rRaw,BS,Nch,[]);
			rRaw = permute(rRaw,[1 3 2]);
			rRaw = reshape(rRaw,[],Nch);
		end
	else
		rRaw = dlmread([src_dir filesep Fnm],delimI); %reads data file
		if TS	%If timestamp is present in files
			rRaw = rRaw(:,2:end);	%removes timestamp
		end
		S = length(rRaw(:,1));
		if round(S/BS)~=S/BS
			error(['Data is not an integer number of blocks'])
		end
	end
    [f,SPL] = calcSPL_v2(rRaw,sampleRate,BS,WNDO); clear rRaw f;
    
    dB = ref+ 10*log10(SPL*ChRef); clear SPL;   %Converts into dB using calibration and reference values
    dB = dB(DPoints:end-DPoints,:);   %removes junk points at beginning and end
     
    if ~strcmpi('None',atmAbs)
        %Corrects for absorption by atmosphere - Spectra are now representative 
        %of sound propagated to microphone distance with no absorption.
        [alpha,alpha_s] = atmAbsorption_v2(fx,AT,AP,RH,atmAbs); 
        [x,y] = meshgrid(ChD,alpha); dB = dB + x.*y;    %Removes absorption at current conditions
        clear alpha x y;
    end
    
    if strcmpi('Radial',NDs)    %Correct to normalized distance using simple radial propagation
        if sum(abs(MO))~=0
            dBCorrected = micOrientationCorrect(fx,MO,dB); %Corrects spectra for microphone orientation
        else
            dBCorrected = dB;
        end
        [dBCorrected, al] = NormDistance(dBCorrected,ChD/D,NormD,dBShift);
        if strcmpi('Standard',atmAbs)   %Adds absorption at standard day conditions if selected
            [x,y] = meshgrid(ChD.*al',alpha_s); dBCorrected = dBCorrected - x.*y;
        end
    elseif strcmpi('NS',NDs)        %Correct to normalized distance using noise source location

        [dBCorrected, al, ph] = NormDistance(dB,Std,ChD/D,ChPol*pi/180,NormD,dBShift);
        if sum(abs(MO))~=0
            MOph = abs(ph*180/pi - repmat(MO+ChPol,length(Std),1));
            dBCorrected = micOrientationCorrect(repmat(fx,1,length(MO)),MOph,dBCorrected); %Corrects spectra for microphone orientation        
        end
        if strcmpi('Standard',atmAbs)
            [x,y] = meshgrid(ChD,alpha_s); dBCorrected = dBCorrected - x.*y.*al;
        end
    else
        if sum(abs(MO))~=0
            dBCorrected = micOrientationCorrect(fx,MO,dB);
        else
            dBCorrected = dB;
        end
        if strcmpi('Standard',atmAbs)   %Adds absorption at standard day conditions if selected
            [x,y] = meshgrid(ChD,alpha_s); dBCorrected = dBCorrected - x.*y;
        end
    end
    clear al x y dB alpha_s MOph;
    
    dBCorrected_DT = dBCorrected;
    if DT~=0
        for m = 1:Nch   %Smoothes and removes forcing harmonics from spectrum. If FH==0 (forcing frequency) program doesn't look for harmonics
            dBCorrected_DT(:,m) = ToneSubtract_v5(fx,dBCorrected(:,m),'Poly',30);
        %If tones are broadband and poly fitting is failing, use this alternate call:
%             dBCorrected_DT(:,m) = ToneSubtract_v5(fx,dBCorrected(:,m),'Moving',FH,10);
        end
    else
        for m = 1:Nch   %Smoothes forcing harmonics from spectrum.
            dBCorrected_DT(:,m) = ToneSubtract_v5(fx,dBCorrected(:,m),'Poly',30);
        %If tones are broadband and poly fitting is failing, use this alternate call:
%             dBCorrected_DT(:,m) = ToneSubtract_v5(fx,dBCorrected(:,m),'Moving',0,10);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OUTPUT FILE GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%  PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure  %Plots all spectra on a single figure
	CM = colormap('jet');
    semilogx(Std,dBCorrected(:,1),'Color',CM(round(64/Nch*1),:));
    LBL = 'Ch 1';  %Variable contains labels which will become column headers on output files
    hold on
    for m = 2:Nch
        semilogx(Std,dBCorrected(:,m),'Color',CM(round(64/Nch*m),:));
        LBL = [LBL '\tCh ' num2str(m)]; %#ok<AGROW>
    end
    LEGN = mat2cell(ChPol,1,ones(1,Nch));
    for m = 1:Nch   %Plots smoothed spectra on top of original spectra
        semilogx(Std,dBCorrected_DT(:,m),'--k');
        LEGN{m} = num2str(LEGN{m});
    end
    legend(LEGN,'Location','NorthWest'); clear LEGN;
    xlabel(['Strouhal Number for D = ' num2str(D) ' m Jet at M = ' num2str(M)]);
    ylabel('SPL (dB)');
    title({['Average Spectrum for: ' Fnm],['Data Normalized for ' num2str(NormD) ' x/D'],['Reynolds Number: ' num2str(ReN,'%.0f')]});
    grid on;
    if length(fLimits)==2
        xlim(fLimits);
    else
        axis(fLimits);
    end
    saveas(gcf,[out_dir filesep Fnm(1:end-4) '.fig']);   %Saves figure as .fig
    saveFigure_v2(gcf,[out_dir filesep Fnm(1:end-4)],300)   %Saves figure as .png
    close    %Closes figure
    
    %%%%%  SPECTRAL DATA FILES  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([out_dir filesep Fnm(1:end-3) 'fft' EXT],'w');   %Opens .fft* file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number' delimO LBL '\n']);
    fclose(fid);
    dlmwrite([out_dir filesep Fnm(1:end-3) 'fft' EXT],[fx Std dBCorrected], 'delimiter', delimO, '-append');   %writes spectral data to file
    
    fid = fopen([out_dir filesep Fnm(1:end-3) 'S.fft' EXT],'w');   %Opens .S.fft* file to write column headers
    fprintf(fid,['Frequency (Hz)\tStrouhal Number' delimO LBL '\n']);
    fclose(fid); clear LBL;
    dlmwrite([out_dir filesep Fnm(1:end-3) 'S.fft' EXT],[fx Std dBCorrected_DT], 'delimiter', delimO, '-append');   %writes smoothed spectral data to file

    %%%%%  SET FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([out_dir filesep Fnm(1:end-3) 'SET'],'w');
    fprintf(fid,'PROCESS VALUES FOR DATA\n');
    fprintf(fid,['Processing Time Stamp:' delimO datestr(now,0) '\n']);
    fprintf(fid,['Sampling Rate (Hz):' delimO num2str(sampleRate) '\n']);
    fprintf(fid,['Block Size:' delimO num2str(BS) '\n']);
    fprintf(fid,['Reference Amplitude (dB):' delimO num2str(ref) '\n']);
    fprintf(fid,['Window:' delimO WNDO '\n']);
    fprintf(fid,['Normalization Distance (x/D):' delimO num2str(NormD) '\n']);    
    fprintf(fid,['Normalization Method:' delimO NDs '\n']);
    fprintf(fid,['Number of Channels:' delimO num2str(Nch) '\n']);
        %Records location and name of each calibration file used.
    fprintf(fid,['Location of Calibration Files:' delimO CalPathP '\n']);
    fprintf(fid,'Calibration Files:');
    for nn = 1:Nch
        fprintf(fid,[delimO CalFiles{nn}]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mic Distance (m):'); fprintf(fid,'\t%.2f',ChD); fprintf(fid,'\n');
    fprintf(fid,'Mic Polar Angle (deg):'); fprintf(fid,'\t%.2f',ChPol); fprintf(fid,'\n');
    fprintf(fid,'Mic Orientation Angle (deg):'); fprintf(fid,'\t%.2f',MO); fprintf(fid,'\n');
    if DT==0
        fprintf(fid,'Tone Subtraction:\tNO\n');
    else
        fprintf(fid,'Tone Subtraction:\tYES\n');
    end
    fprintf(fid,['Spectrum offset (dB):' delimO num2str(dBShift) '\n']);
    fprintf(fid,['Trim Ends of Spectrum (N):' delimO num2str(DPoints) '\n']);
    if strcmpi('Standard',atmAbs)
        fprintf(fid,'Atmospheric Absorption:\tStandard Day\n');
    elseif strcmpi('Lossless',atmAbs)
        fprintf(fid,'Atmospheric Absorption:\tLossless Propagation\n');
    else
        fprintf(fid,'Atmospheric Absorption:\tNo Correction\n');
    end
    fprintf(fid,['Jet Diameter (m):' delimO num2str(D) '\n']);
    fprintf(fid,['Mach Number:' delimO num2str(M) '\n']);
    fprintf(fid,['Ambient Temperature (K):' delimO num2str(AT+273.15) '\n']);
    fprintf(fid,['Ambient Pressure (kPa):' delimO num2str(AP) '\n']);
    fprintf(fid,['Relative Humidity (%%):' delimO num2str(RH) '\n']);
    fprintf(fid,['Stagnation Temperature (K):' delimO num2str(To) '\n']);
    fprintf(fid,['Exit Temperature (K):' delimO num2str(Temp) '\n']);
    fprintf(fid,['Total Temperature Ratio:' delimO num2str(To/(AT+273.15)) '\n']);
    fprintf(fid,['Exit Velocity (m/s):' delimO num2str(M*sqrt(1.4*287.05*Temp)) '\n']);
    fprintf(fid,['Jet Velocity Ratio (Uj/a_inf):' delimO num2str(M*sqrt(Temp/(AT+273.15))) '\n']);
    fprintf(fid,['Exit Density (kg/m^3):' delimO num2str(1.01e5/(287.05*Temp)) '\n']);
    fprintf(fid,['Exit Viscosity (Pa s):' delimO num2str((1.827e-5*(0.555*291.15+66.667)/(0.555*Temp+66.667).*(Temp/291.15).^(3/2))) '\n']);   %equation is Sutherlands formula - 291.15 is reference temperature
    fprintf(fid,['Jet Reynolds Number:' delimO num2str(ReN,'%.0f') '\n']);
    
        %Calculates and prints OASPL for all channels    
    Bg = find(Std > iLimits(1),1,'first'); %Excludes spectral components below 0.04 Strouhal number
    Ed = find(Std > iLimits(2),1,'first');    %Excludes spectral components above 4 Strouhal number
    if isempty(Bg)
        Bg = 1;
    end
    if isempty(Ed)
        Ed = length(Std);
    end
    OASPL = 10*log10(trapz(fx(Bg:Ed),10.^(dBCorrected(Bg:Ed,:)/10),1)); clear dBCorrected;
    fprintf(fid,'OASPL (dB):'); fprintf(fid,'\t%.3f',OASPL); fprintf(fid,'\n');
        %Calculates and prints the average acoustic energy
    delta_Pol = abs(diff(ChPol)); IW = ([delta_Pol/2 0] + [0 delta_Pol/2])/(max(ChPol)-min(ChPol));
    AAE = 10*log10(sum(IW.*10.^(OASPL/10)));
    fprintf(fid,['Average Energy (dB):' delimO num2str(AAE,'%.3f') '\n']);
        %Calculates and prints the normalized jet total power - assumes
        %azimuthally symmetric and fully expanded. The solid angle of the 
        %integration is divided out giving (W/sr). This power per solid
        %angle is then normalized by the jet mechanical power per solid
        %angle of a complete sphere (i.e. 4pi).
    NJP = 10*log10(32*NormD^2*pref^2*sqrt((AT+273.15)/Temp)*...
        trapz(ChPol*pi/180,10.^(OASPL/10).*sin(ChPol*pi/180))/...
        ((1.4*AP*1000)^2*M^3*(cos(ChPol(1)*pi/180)-cos(ChPol(end)*pi/180))));
    fprintf(fid,['Normalized Jet Power (dB):' delimO num2str(NJP,'%.3f') '\n']);
    
        %Calculates and prints the peak SPL magnitude and Strouhal number
    fprintf(fid,'Peak of SPL Curve\n');    
    fprintf(fid,'  Magnitude (dB):');
    for nn = 1:Nch
        [c,I] = max(dBCorrected_DT(1:Ed,nn));
        I(I < 6) = 6;
        c = mean(dBCorrected_DT(I-5:I+5,nn));
        fprintf(fid,[delimO num2str(c,'%.3f')]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'  Frequency (Std):');
    for nn = 1:Nch
        [c,I] = max(dBCorrected_DT(1:Ed,nn));
        fprintf(fid,[delimO num2str(Std(I),'%.4f')]);
    end
    fprintf(fid,'\n');
        %Calculates and prints the high frequency slope of the spectra 
        %including uncertainty as seen on semilogx scale.
    fprintf(fid,'Slope of SPL Curve (Strouhal axis):');
    for nn = 1:Nch
        I = find(Std >= 0.5,1,'first');   %Ignores all points to left of Std = 0.7
        Pval = mmpolyfit(log10(Std(I:Ed)),dBCorrected_DT(I:Ed,nn),1,'Weight',1./Std(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
            %Calculates uncertainty of fit
        Serr = sum((Pval(1)*log10(Std(I:Ed))+Pval(2)-dBCorrected_DT(I:Ed,nn)).^2);
        merr = length(Std(I:Ed));
        Derr = merr*sum(log10(Std(I:Ed)).^2) - sum(log10(Std(I:Ed)))^2;
        PE = sqrt(Serr*merr/(merr-2)/Derr);
        
        fprintf(fid,[delimO num2str(Pval(1),'%10.3f') ' +/- ' num2str(PE,'%10.3f')]);
    end
    fprintf(fid,'\n');
    
        %If forcing is present, calculates/prints several relevant values
    if SHarmonics
        fprintf(fid,['Forcing Frequency (Hz):' delimO num2str(FH) '\n']);
        fprintf(fid,['Forcing Strouhal Number:' delimO num2str(StDF) '\n']);
            %Calculates and prints OASPL for all channels
        if DT~=0
            OASPL_DT = 10*log10(trapz(fx(Bg:Ed),10.^(dBCorrected_DT(Bg:Ed,:)/10),1)); clear dBCorrected_DT;
            fprintf(fid,'OASPL_Detoned (dB):'); fprintf(fid,'\t%.3f',OASPL_DT); fprintf(fid,'\n');
                %Calculates and prints the average acoustic energy
            AAE_DT = 10*log10(sum(IW.*10.^(OASPL_DT/10)));
            fprintf(fid,['Average Energy_Detoned (dB):' delimO num2str(AAE_DT,'%.3f') '\n']);
                %Calculates and prints the normalized jet total power
            NJP_DT = 10*log10(32*NormD^2*pref^2*sqrt((AT+273.15)/Temp)*...
                trapz(ChPol*pi/180,10.^(OASPL_DT/10).*sin(ChPol*pi/180))/...
                ((1.4*AP*1000)^2*M^3*(cos(ChPol(1)*pi/180)-cos(ChPol(end)*pi/180))));
            fprintf(fid,['Normalized Jet Power_Detoned (dB):' delimO num2str(NJP_DT,'%.3f') '\n']);     
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
                
                if length(BT) > 1
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
            end
            
                %Calculates and prints delta OASPL
            fprintf(fid,'dOASPL (dB):'); fprintf(fid,'\t%.4f',OASPL - OASPLB); fprintf(fid,'\n');
                %Calculates and prints change in average acoustic energy
            fprintf(fid,['dAAE (dB):' delimO num2str(AAE-AAEB,'%.4f') '\n']);
                %Calculates and prints change in normalized jet power
            fprintf(fid,['dNJP (dB):' delimO num2str(NJP-NJPB,'%.4f') '\n']);
            
            if DT~=0
                    %Calculates and prints delta OASPL
                fprintf(fid,'dOASPL_Detoned (dB):'); fprintf(fid,'\t%.4f',OASPL_DT - OASPLB_DT); fprintf(fid,'\n');
                    %Calculates and prints change in average acoustic energy
                fprintf(fid,['dAAE_Detoned (dB):' delimO num2str(AAE_DT-AAEB_DT,'%.4f') '\n']);
                    %Calculates and prints change in normalized jet power
                fprintf(fid,['dNJP_Detoned (dB):' delimO num2str(NJP_DT-NJPB_DT,'%.4f') '\n']);
            end
        end       
    else
            %Calculates and prints OASPL for all channels    
        OASPL_DT = 10*log10(trapz(fx(Bg:Ed),10.^(dBCorrected_DT(Bg:Ed,:)/10),1)); clear dBCorrected;
        fprintf(fid,'OASPL_Smoothed (dB):'); fprintf(fid,'\t%.3f',OASPL_DT); fprintf(fid,'\n');
            %Calculates and prints the average acoustic energy
        AAE_DT = 10*log10(sum(IW.*10.^(OASPL_DT/10)));
        fprintf(fid,['Average Energy_Smoothed (dB):' delimO num2str(AAE_DT,'%.3f') '\n']);
            %Calculates and prints the normalized jet total power
        NJP_DT = 10*log10(32*NormD^2*pref^2*sqrt((AT+273.15)/Temp)*...
            trapz(ChPol*pi/180,10.^(OASPL_DT/10).*sin(ChPol*pi/180))/...
            ((1.4*AP*1000)^2*M^3*(cos(ChPol(1)*pi/180)-cos(ChPol(end)*pi/180))));
        fprintf(fid,['Normalized Jet Power_Smoothed (dB):' delimO num2str(NJP_DT,'%.3f') '\n']);
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECORD LOG INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('..Recording Log Information')
dep_dir = [out_dir filesep 'Dependencies'];
mkdir(dep_dir)

copyfile(PARAMS,[dep_dir filesep 'PARAMS.mat'])

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

