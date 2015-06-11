function flist = getFList(path,ext,subsetFlag,subFlagString)
% flist = getFList(path);
% flist = getFList(path,ext);
% flist = getFList(path,ext,1);
% flist = getFList(path,ext,2,subFlagString);
%
%
% CALLS
% 1) flist = getFList(path); - returns the complete list of files and dirs 
%		in the location "path".
% 2) flist = getFList(path,ext); - returns the set of files with extension
%		"ext" unless "ext"=='dir' in which case a list of directories is
%		returned.
% 3) flist = getFList(path,ext,1); - same as second call, but presents the
%		user with a dialog box to select a subset of the first pass result.
% 4) flist = getFList(path,ext,2,subFlagString); - same as second call, but
%		uses the string in variable "subFlagString" to automatically select
%		a subset of the first pass result.
%
%
% path = 'C:\whatever';
% ext = '.junk'; - selects for directories if ext=='dir'
% subsetFlag = 1; %1 = manual subset, 2 = auto subset from argument




flist = struct2cell(dir(path));
flist = flist(:,3:end);	%Eliminates directory back links '.' & '..'

if exist('ext','var')
	if strcmpi(ext(1),'.') %Strips period from extension if it exists
		ext = ext(2:end);
	end
	
	if strcmpi(ext,'dir')
		%Selects only directories
		flist = flist(1,cell2mat(flist(4,:)));
	else
		%Eliminates all files not having extension "ext"
		flist = flist(1,:);
		keep = zeros(size(flist));
		for n = 1:length(flist)
			keep(n) = strcmpi(flist{n}(end-length(ext)+1:end),ext);
		end
		flist = flist(logical(keep));
	end

	if exist('subsetFlag','var')
		switch subsetFlag
			case 1
					%Lets user choose the subset of files in that directory for gathering.
				[keep,ok] = listdlg('PromptString','Select Desired Files',...
								'SelectionMode','multiple',...
								'ListString',flist);
				if ok==0
					error('getFList Terminated Due to User Selection of Cancel')
				end
			case 2
					%Uses string specified in call to autoselect a subset.
				if exist('subFlagString','var')
					q = strfind(flist,subFlagString);
					keep = ~cellfun('isempty',q);
				else
					error('subFlagString must exist')
				end
			otherwise
				keep = true(size(flist));
		end
		flist = flist(keep);
	end
else
	flist = flist(1,:);
end