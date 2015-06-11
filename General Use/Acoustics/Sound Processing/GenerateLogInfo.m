function GenerateLogInfo(functionname,tstamp,PARAMS,pp,dep_dir)           
    disp('..Recording Log Information');
    copyfile(PARAMS,[dep_dir '\' 'PARAMS.mat']);
    for n = 1:pp.Nch   %Copies all calibration files to dependencies folder
        copyfile([pp.CalPath '\' pp.CalFiles{n}],[dep_dir '\' pp.CalFiles{n}])
    end 
    dep = getFileDependencies(functionname);
    for n = 1:length(dep)   %Copies all functions used to dependencies folder
        copyfile(dep{n},[dep_dir '\']);
    end   
    fid = fopen([dep_dir '\' 'ProcessedFiles.txt'],'w');  %Writes log file containing list of files processed 
    fprintf(fid,['Files processed on: ' datestr(tstamp,1) '\n']);
    fprintf(fid,['Start Time: ' datestr(tstamp,13) ', End Time: ' datestr(now,13) '\n']);
    for n = 1:length(pp.flist)
        fprintf(fid,[pp.flist{n} '\n']);
    end
    fclose(fid);
end