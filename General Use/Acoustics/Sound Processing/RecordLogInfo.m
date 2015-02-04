function RecordLogInfo(functionname,tstamp,PARAMS,pp,dep_dir,method)
%Records log information.  Includes the processing parameters used, a list
%of files that were processed, a copy of the calibration files that were
%used, and a copy of all the .m files that were used.  Files are copied
%into the directory specified by dep_dir.  If the user desires, files can
%be compressed using the zip method by adding either '-c' or '-compress'
%for the input variable 'method', otherwise this variable is unnecessary.

%Last updated by Michael Crawley on 2011-07-04

    disp('..Recording Log Information');
    if ~exist('method','var')
       method = []; 
    end
    dep = getFileDependencies(functionname);
     
    if strcmpi(method,'-c') || strcmpi(method,'-compress') %compress files using zip format    
        fid = fopen('ProcessedFiles.txt','w');  %Writes temporary log file containing list of files processed 
        fprintf(fid,['Files processed on: ' datestr(tstamp,1) '\n']);
        fprintf(fid,['Start Time: ' datestr(tstamp,13) ', End Time: ' datestr(now,13) '\n']);
        for n = 1:length(pp.flist)
            fprintf(fid,[pp.flist{n} '\n']);
        end
        fclose(fid);        
        CalFiles = cell(1,pp.Nch);
        for n = 1:pp.Nch
            CalFiles{n} = fullfile(pp.CalPath,pp.CalFiles{n}); %adds full pathname to calfiles
        end        
        zip(dep_dir,[dep', PARAMS, CalFiles, 'ProcessedFiles.txt']); %compresses files into specified folder
        delete('ProcessedFiles.txt'); %deletes temp file
    else %standard copy
        mkdir(dep_dir);        
        copyfile(PARAMS,[dep_dir '\' 'PARAMS.mat']);
        for n = 1:pp.Nch   %Copies all calibration files to dependencies folder
            copyfile([pp.CalPath '\' pp.CalFiles{n}],[dep_dir '\' pp.CalFiles{n}])
        end
    
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
end