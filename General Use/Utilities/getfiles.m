function [flist,src_dir,folders] = getfiles(condition,varargin)
%[flist,src_dir,folders] = getfiles(condition,varargin)
%Function returns list of files and folders (and source directory) based on
%a given condition.
%Inputs:
%           condition:  conditional string on which returned files are 
%                       dependent.  Accepts wildcards.  [] returns all 
%                       files in directory.
%           options:
%                       path:   specify top-level directory for searching.
%                               If this is not specified, the program 
%                               prompts the user.
%                       '-a' or '-all': searches subdirectories in addition
%                                       to top-level.
%                       '-s' or '-sub': allows user to select a subset of
%                                       the found files.

    if ~exist('condition','var'), condition = ''; end 
    if isempty(varargin), varargin = {''}; end
    
    %Determine search directory
    test = cellfun(@(lst) isdir(lst),varargin);
    if any(test)
        src_dir = varargin{test};
    else
        src_dir = uigetdir('Specify Directory'); 
    end
    
    %find files within directory that match given condition
    flist = findfiles(src_dir,condition);
    
    %find folders and contents within
    if any(strcmpi(varargin,'-all')) || any(strcmpi(varargin,'-a'))  
        tmp = genpath(src_dir);
        cutoff = length(src_dir)+2; %cutoff for folder name
        folders = regexp(tmp,':','split'); %find subdirs
        folders = folders(2:end); %remove top level directory from subdir list
        folder_keep = true(1,length(folders)); %logical check for folder contents
        for n = 1:length(folders)
           tflist = findfiles(folders{n},condition); %find files in given subdir
           if ~isempty(tflist)
               tflist = cellfun(@(file) fullfile(folders{n}(cutoff:end),file),tflist,'UniformOutput',false); %append subdir name to file
               flist = [flist, tflist]; %append new subdir filelist to dir filelist
           else 
               folder_keep(n) = false;
           end
        end
        folders = folders(folder_keep);
    else
        folders = cell(0);
    end
    
    %Lets user choose the subset of files in that directory for gathering
    if any(strcmpi(varargin,'-sub')) || any(strcmpi(varargin,'-s')) 
        l = max(cellfun(@length,flist));
        [keep,ok] = listdlg('PromptString','Select Desired Files','SelectionMode','multiple','ListSize',[8*l 300],'ListString',flist);
        if ok==0
            error('Function Terminated Due to User Selection of Cancel')
        end
        flist = flist(keep);
    end
    flist = flist';
end

function [tflist] = findfiles(src_dir,condition)
    list = dir([src_dir,filesep,condition]);
    chk = [list.isdir];
    tflist = {list(~chk).name};
end