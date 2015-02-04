function [flist,src_dir,folders] = getfiles(condition,varargin)
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
    list = dir([src_dir,filesep,condition]);
    chk = [list.isdir];
    flist = {list(~chk).name};
    
    %find folders and contents within
    tmp = dir(src_dir);
    chk = [tmp.isdir];       
    folders = {tmp(chk).name};
    folders = folders(3:end);
    if any(strcmpi(varargin,'-all')) || any(strcmpi(varargin,'-a'))  
       toempty = false(1,length(folders));
       sflist = {}; sfolders = {};
       for n = 1:length(folders)
          [tflist,t_dir,tfolders] = getfiles(condition,[src_dir filesep folders{n}],'-a');
          tflist = cellfun(@(file) fullfile(folders{n},file),tflist,'UniformOutput',false);
          tfolders = cellfun(@(folder) fullfile(t_dir,folder),tfolders,'UniformOutput',false);
          sflist = [sflist,tflist]; %append new files
          if ~isempty(tfolders)
              sfolders = [sfolders,tfolders];
          end
          if isempty(tfolders) && isempty(tflist)
              toempty(n) = true;            
          end
       end
       folders = folders(~toempty);
       flist = [flist,sflist];
       folders = cellfun(@(folder) fullfile(src_dir,folder),folders,'UniformOutput',false);
       folders = [folders,sfolders];
    else
        folders = cellfun(@(folder) fullfile(src_dir,folder),folders,'UniformOutput',false);
    end
    
    %Lets user choose the subset of files in that directory for gathering
    if any(strcmpi(varargin,'-sub')) || any(strcmpi(varargin,'-s'))        
        [keep,ok] = listdlg('PromptString','Select Desired Files',...
                        'SelectionMode','multiple',...
                        'ListString',flist);
        if ok==0
            error('getFList Terminated Due to User Selection of Cancel')
        end
        flist = flist(keep);
    end    
end

%there is an error in this code; the folder names are incorrect