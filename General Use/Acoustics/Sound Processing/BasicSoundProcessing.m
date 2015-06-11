function BasicSoundProcessing(PARAMS)
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

%Last updated by Michael Crawley on 2011-07-21

%INPUTS
%  PARAMS - The path to the .mat file containing the parameters. This file 
%  must always be used, but some of the variables in it are optional. See
%  GenerateProcessingParams.m for more information regarding the
%  requirements for the PARAMS .mat file.

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


    %%% Load Processing Parameters
    pp = load(PARAMS);
    pp = GenerateProcessingParams(pp);
    pp = ProcessCalFiles(pp);
    
    %%% Create save folders 
    tstamp = now;
    out_dir = [pp.src_dir '\' datestr(tstamp,30)];   
    dep_dir = [out_dir '\' 'Dependencies'];
    mkdir(out_dir);
    
    %%% Process baseline Files
    for n = 1:pp.nBaselines
        filename = pp.flist{n};
        data = ReadFileName(filename,pp);
        rRaw = ReadFile(filename,pp);
        data = CalcSpectra(rRaw,pp,data);
        results = CalcResults(data,pp);
        GenerateFFTFile(filename,data,pp,out_dir);
        GenerateSETFile(filename,results,data,pp,out_dir);
        GenerateFigureFile(filename,data,pp,out_dir);
    end
    
    %%% Pull baseline results
    if pp.nBaselines ~= length(pp.flist) && pp.nBaselines ~= 0
        pp = ProcessBaselineResults(pp,out_dir);
    end

    %%% Process forced Files
    for n = n+1:length(pp.flist)
        filename = pp.flist{n};
        data = ReadFileName(filename,pp);
        rRaw = ReadFile(filename,pp);
        data = CalcSpectra(rRaw,pp,data);
        results = CalcResults(data,pp);
        GenerateFFTFile(filename,data,pp,out_dir);
        GenerateSETFile(filename,results,data,pp,out_dir);
        GenerateFigureFile(filename,data,pp,out_dir);
    end
    
    %%% Record Log Information
    RecordLogInfo(mfilename,tstamp,PARAMS,pp,dep_dir);    
end