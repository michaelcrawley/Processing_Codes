function [pp] = GenerateProcessingParams(pp,suppress) 
%Checks processing parameter structure (pp) for required fields, and asks
%user for input for any missing fields.  Updates processing parameter
%structures.  Miscellaneous fields are added to the structure which are not
%expected to change between processing data sets.

%Last updated by Michael Crawley on 2011-09-09

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

    pp.TS = false;		%false == no timestamp present in data files
    pp.TS_cal = false;	%false == no timestamp present in calibration files
    pp.delimI = '\t';	%Delimiter to check for in file reading 
    pp.delimO = '\t';	%Delimiter for file writing - As of 20101006 other functions assume output files are tab delimited...change at own risk.
    pp.fLimits = [0.01 8]; %The limits of the plotting on the figure windows. Use [xx xx] to limit only abscissa. Use [xx xx xx xx] to limit all dimensions.
    pp.iLimits = [0.04 4]; %The limits of integration (in Strouhal number) for determining OASPL.
    pp.pref = 20e-6;   %Pascals - reference pressure - used in calculating normalized jet power

    if strcmpi('None',pp.NDs)  %Forces Normalization distance to invalid if no propagation correction is wanted
        pp.NormD = NaN;
    end
    if isfield(pp,'src_dir') && isfield(pp,'flist') && isfield(pp,'CalPath') && isfield(pp,'CalFiles')
        pp.SD = 'N.A.';
    end
    
    if ~exist('suppress','var')
        suppress = 0;
    elseif strcmpi(suppress,'-suppress') || strcmpi(suppress,'-s') || strcmpi(suppress,'-noquestion')
        suppress = 1;
    else
        suppress = 0;
    end
    

        %Asks user to verify proper parameters
    disp('Current Processing Values: ');
    disp(['                 Ambient Temp (C): ' num2str(pp.AT)]);
    disp(['           Ambient Pressure (kPa): ' num2str(pp.AP)]);
    disp(['            Relative Humidity (%): ' num2str(pp.RH)]);
    disp(['                 Jet Diameter (m): ' num2str(pp.D)]);
    disp(['                      Mach Number: ' num2str(pp.M)]);
    disp(['      Reference Signal Level (dB): ' num2str(pp.ref)]);
    disp(['               Sampling Rate (Hz): ' num2str(pp.sampleRate)]);
    disp(['                Samples Per Block: ' num2str(pp.BS)]);
    disp(['Mic Normalization Distance Method: ' pp.NDs]);
    disp([' Mic Normalization Distance (r/D): ' num2str(pp.NormD)]);
    disp(['       Mic Polar Angles (degrees): ' num2str(pp.ChPol,'%.0f ')]);
    disp(['  Mic Orientation Angle (degrees): ' num2str(pp.MO,'%.0f ')]);
    disp(['                 Mic Distance (m): ' num2str(pp.ChD,'%.2f ')]);
    if pp.DT==0
        disp('                 Tone Subtraction: NO');
    else 
        disp('                 Tone Subtraction: YES');
    end
    disp(['    Spectrum Vertical Offset (dB): ' num2str(pp.dBShift)]);
    disp(['       Trim ends of sprectrum (N): ' num2str(pp.DPoints)]);
    disp(['                      Data Window: ' num2str(pp.WNDO)]);
    disp(['                   Atm Absorption: ' pp.atmAbs]);
    disp(['               Starting Directory: ' pp.SD]);
    if isfield(pp,'src_dir') 
        disp(['            Data Source Directory: ' pp.src_dir]); 
    end
    if isfield(pp,'flist')
        disp('                   Data File List: PRESENT');
    else
        disp('                   Data File List: NOT PRESENT');
    end
    if isfield(pp,'CalPath')
        disp(['       Calibration File Directory: ' pp.CalPath]); 
    end
    if isfield(pp,'CalFiles')
        disp('            Calibration File List: PRESENT');
    else
        disp('            Calibration File List: NOT PRESENT');
    end
    
    if ~suppress
        yn = questdlg_timeout('Is the information correct? - Program assumes yes if no user selection after 60 seconds.','Parameter Confirmation',60);
        if ~strcmpi(yn,'Yes') && ~strcmpi(yn,'')
            error('Fix Values in code and run again')
        end
        if strcmpi(yn,'')
            disp('  NO user confirmation of parameters...Assuming full-auto operation')
        end
    end

    if ~isfield(pp,'src_dir')
        pp.src_dir = uigetdir(pp.SD,'Specify Directory Containing Data'); 
    end
    USelect = false;    %switch variable that allows user to select a file list subset if the file list was not submitted in the "PARAMS" file
    if ~isfield(pp,'flist')
        USelect = true;
        flist = struct2cell(dir(pp.src_dir));
        tmp = strfind(flist(1,:),['.' pp.EXT]); %extracts list of data files - ignores all other files
        keep = ~cellfun('isempty',tmp); clear tmp;
        pp.flist = flist(1,logical(keep));
    end
    
    %make sure flist is a cell (necessary if only one file is being
    %processed)
    if ~iscell(pp.flist)
       pp.flist = cellstr(pp.flist); 
    end
    
    pp.flist = IDBaselines(pp.src_dir,pp.flist);    %rewrites file names of Baselines
    if USelect
        [keep,ok] = listdlg('PromptString','Select files to process:',...
                        'SelectionMode','multiple',...
                        'ListString',pp.flist);    %Asks the user to select the subset of data files they wish to process
        if ok==0
            error('Program Terminated Due to User Selection of Cancel')
        end
        pp.flist = pp.flist(keep); clear keep ok;
    end
    clear USelect;

    if ~isfield(pp,'CalPath') || ~isfield(pp,'CalFiles')
        %Asks user to select the calibration files to use and verifies proper
        %naming convention.
        [pp.CalFiles,pp.CalPath] = uigetfile(['*.' pp.EXT],'Select Calibration File(s) for Data',pp.SD,'MultiSelect','on');
    end
    
    %make sure CalFiles is a cell
    if ~iscell(pp.CalFiles)
        pp.CalFiles = cellstr(pp.CalFiles);
    end
    
    pp.Nch = length(pp.CalFiles); %Determines the number of channels of data by the number of calibration files
    if length(pp.ChPol)~=pp.Nch
        error('    Microphone information mismatch!! Number of calibration files does not match microphone location information.')
    end
    
    pp.fx = (0:ceil((pp.BS+1)/2)-1)'*(pp.sampleRate/pp.BS);   %Frequency axis
    pp.fx = pp.fx(pp.DPoints:end-pp.DPoints); %First and last points will be thrown out as garbage    
    
    %Reorders file list so that baseline files are before forced
    %files.  Requires that the baseline files have the name 'Baseline' in them. 
    match = strfind(pp.flist,'Baseline');
    forced = cellfun(@isempty,match);
    pp.flist = [pp.flist(~forced)  pp.flist(forced)];
    pp.nBaselines = length(pp.flist) - sum(forced);
    
    %determines if data files are ascii or binary
    pp.isbinary = true;
    fid = fopen([pp.src_dir '\' pp.flist{1}],'r');
    rRaw = fread(fid,500,'*char'); fclose(fid);	%Reads first 200 bytes assuming they are ascii
    if strcmp('\t',pp.delimI)
        if length(strfind(rRaw',char(9))) > 5	%If delimiter is present, file is not binary
            pp.isbinary = false;
        end
    else
        if length(strfind(rRaw',pp.delimI)) > 5
            pp.isbinary = false;
        end
    end
    if pp.isbinary
        if length(strfind(rRaw',char(10))) > 5
            pp.isbinary = false;
        end
    end
    
    if ~isfield(pp,'BT') %check to see if broadband tones switch is present, if not, set switch to false
       pp.BT = 0; 
    end
end