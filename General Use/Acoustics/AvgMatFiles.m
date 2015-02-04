function [] = AvgMatFiles(flist,endfilename)
%function averages listed matfile data from soundprocessing code
%file names listed in flist must be full file names, or the file must be in
%the working directory or path.  Assumes all files are for forced cases,
%and therefore have 'delta' metrics.
%Currently under construction.  Resultant output matfile is compatible with
%AcousticPlotter for -doaspl and -NLM commands.


    for n = 1:length(flist)
       tmp = load(flist{n});
       [~,filename] = fileparts(flist{n});
       filename = filename(1:end-4); %temp fix for misnamed .mat files missing azimuthal angle information
       Stdf(n) = eval(['tmp.',filename,'.data.StDF']);
       ChPol(:,n) = eval(['tmp.',filename,'.processing_params.ChPol']);
       dOASPL(:,:,n) = eval(['tmp.',filename,'.results.dOASPL']);
       dOASPL_DT(:,:,n) = eval(['tmp.',filename,'.results.dOASPL_DT']);
       dAAE_DT(:,:,n) = eval(['tmp.',filename,'.results.dAAE_DT']);
       Pskewness(:,:,n) = eval(['tmp.',filename,'.results.Pskewness']);
       dPskewness(:,:,n) = eval(['tmp.',filename,'.results.dPskewness']);
       Pkurtosis(:,:,n) = eval(['tmp.',filename,'.results.Pkurtosis']);
       dPkurtosis(:,:,n) = eval(['tmp.',filename,'.results.dPkurtosis']);
    end
    
    %calc averages of metrics
    eval([endfilename,'.data.StDF = mean(Stdf);']);
    eval([endfilename,'.processing_params.ChPol = mean(ChPol,2);']);
    eval([endfilename,'.results.dOASPL = mean(dOASPL,3);']);
    eval([endfilename,'.results.dOASPL_DT = mean(dOASPL_DT,3);']);
    eval([endfilename,'.results.dAAE_DT = mean(dAAE_DT,3);']);
    eval([endfilename,'.results.Pskewness = mean(Pskewness,3);']);
    eval([endfilename,'.results.dPskewness = mean(dPskewness,3);']);
    eval([endfilename,'.results.Pkurtosis = mean(Pkurtosis,3);']);
    eval([endfilename,'.results.dPkurtosis = mean(dPkurtosis,3);']);
    
    %calc std of metrics
    eval([endfilename,'.data.StDFstd = std(Stdf);']);
    eval([endfilename,'.processing_params.ChPolstd = std(ChPol,1,2);']);
    eval([endfilename,'.results.dOASPLstd = std(dOASPL,1,3);']);
    eval([endfilename,'.results.dOASPL_DTstd = std(dOASPL_DT,1,3);']);
    eval([endfilename,'.results.dAAE_DTstd = std(dAAE_DT,1,3);']);
    eval([endfilename,'.results.Pskewnessstd = std(Pskewness,1,3);']);
    eval([endfilename,'.results.dPskewnessstd = std(dPskewness,1,3);']);
    eval([endfilename,'.results.Pkurtosisstd = std(Pkurtosis,1,3);']);
    eval([endfilename,'.results.dPkurtosisstd = std(dPkurtosis,1,3);']);
    
    %add in filenames of files included in average
    eval([endfilename,'.data.flist = flist;']);
    
    %save new matfile
    eval(['save ',endfilename,' ',endfilename]);
end