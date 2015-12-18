function MasterSoundProcessing(PARAMS,cpus)
    %%%Check program inputs
    if ~strcmp(PARAMS(end-3:end),'.mat'), PARAMS = [PARAMS '.mat']; end
    if ~exist('cpus','var'), cpus = FindCoreCount; end    
    
    %%%Open matlab pool & modify loop function, if necessary
    loop = 'for';
    routine = Processing_Routine;
    if cpus > 1
        if matlabpool('size') > 0
            matlabpool close;
        end
        matlabpool(cpus);
        loop = ['par',loop];
    end    
    baserange = ' n = 1:pp.nBaselines, ';
    forcedrange = ' n = pp.nBaselines+1:length(flist), ';    
    
    %%%Generate processing parameters and process calibration files
    pp = GenerateProcessingParams(load(PARAMS),'-s');
    pp = ProcessCalFiles(pp);
    filestamps = gettimestamp(pp.src_dir,pp.flist);
    
    %%% Create save folders 
    tstamp = now;
    if isfield(pp,'save_dir')
        out_dir = [pp.save_dir '\' datestr(tstamp,30)];
    else
        out_dir = [pp.src_dir '\' datestr(tstamp,30)];
    end    
    fig_dir = [out_dir '\Figures'];
    fft_dir = [out_dir '\FFT'];
    dep_dir = [out_dir '\' 'Dependencies'];
    mat_dir = [out_dir '\Mat'];
    set_dir = [out_dir '\SET'];
    mkdir(out_dir);
    mkdir(fig_dir); %create folder for png and fig files
    mkdir(fft_dir); %create folder for fft and S.fft bin files
    mkdir(mat_dir); %create folder for .mat files
    mkdir(set_dir);
    
    %%%Initialize required variables
    flist = pp.flist;    
    
    %%%Process files
    eval([loop,baserange,routine,', end']); %baseline files
    if pp.nBaselines ~= length(flist) && pp.nBaselines ~= 0
        pp = ProcessBaselineResults(pp,set_dir);
    end
    eval([loop,forcedrange,routine,', end']); %forced files
    
    FixFigureFiles(flist,fig_dir);
    
    %%% Record Log Information
    RecordLogInfo(mfilename,tstamp,PARAMS,pp,dep_dir,'-c')
    
    %%% close matlabpool, if necessary
    if cpus > 1
        matlabpool close
    end
    
end

function [routine] = Processing_Routine
    routine = [
        'filename = flist{n};'...
        'data = ReadFileName(filename,pp);'...
        'data.tnum = filestamps(n);'...
        'rRaw = ReadFile(filename,pp);'...
        'data = CalcSpectra(rRaw,pp,data);'...        
        'results = CalcResults(data,pp);'...
        'results = NonlinearMetrics(rRaw,results,pp);'...
        'GenerateFFTFile(filename,data,pp,fft_dir);'...
        'GenerateSETFile(filename,results,data,pp,set_dir);'...
        'GenerateFigureFile(filename,data,pp,fig_dir);'...
        'SaveAcousticData(filename,results,data,pp,mat_dir)'
        ];
end