function [Names,Vals,iI] = gatherSETdata_v2(N,varargin)
%This program allows the user to gather common information from SET files.
%All files must have the same SET file entries since the desired categories
%are only selected once. In order to operate properly, the SET file(s) must 
%be formatted as a tab delimited file in which the first column contains
%the data label for a given row.Upon calling the program, the user is 
%presented with a dialogue box for choosing the directory containing the 
%files. The program then isolates the SET files in that directory and 
%presents them to the user so they may pick the subset they want. Next, the
%program presents the user with a list of the available entries found in 
%the first file to allow the user to select a subset of these. Lastly, if 
%the number of files gathered from the first directory is less than N, the 
%program will present the user with the choices of directory and files 
%until the total count meets/exceeds N.
%
%INPUTS
% N = 21;  %Number of Spectra
% SD = 'My Computer'; %Starting directory
% nFlag = 'm0';	%Auto-select flag
% iI = [1 3 6 7];	%The row numbers of the desired information
%
%OUTPUTS
% Names - a cell array containing the names of the data entries gathered by
% the program.
% Vals - a cell array containing the composited entries gathered by the
% program.
% iI - a vector containing the row numbers of the extracted information.
%
%This function can be called in three forms:
% 1) [Names,Vals,iI] = gatherSETdata(N);  -If no starting directory is given,
% the program uses the default starting point.
% 2) [Names,Vals,iI] = gatherSETdata(N,SD); -the program uses the user
% specified starting directory.
% 2) [Names,Vals,iI] = gatherSETdata(N,SD,nFlag,iI); -the program uses the user
% specified directory and search flag to automatically select the
% set of files from that directory. In this case, the argument "N" is
% ignored.

flag = false;
if isempty(varargin)
    SD = 'My Computer';
elseif nargin==2
    SD = varargin{1};
elseif nargin==4
	SD = varargin{1};
	nFlag = varargin{2}; flag = true;
	iI = varargin{3};
	N = 1;
else
	error('Innappropriate number of arguments');
end

done2 = false;
toto = 0;
SPfile = [];
while ~done2
    
	if flag
		SPpath = SD;
	else
		SPpath = uigetdir(SD,'Specify Directory Containing Data');
	end
	flist = struct2cell(dir(SPpath));
        %Eliminates all files not having .SET in name
    q = strfind(flist(1,:),'.SET');
    keep = ~cellfun('isempty',q);
    flist = flist(1,logical(keep));
	if flag
		q = strfind(flist(1,:),nFlag);
		keep = ~cellfun('isempty',q);
	else
			%Lets user choose the subset of files in that directory for gathering.
		[keep,ok] = listdlg('PromptString','Select files to process:',...
						'SelectionMode','multiple',...
						'ListString',flist);
		if ok==0
			error('Program Terminated Due to User Selection of Cancel')
		end
	end
    SPf = flist(keep);
    
        %Scans through all the kept files extracting the information
    SPfile = [SPfile SPf];
    totn = toto+length(SPf);
	if flag
		Vals = cell(size(iI));
	end
    for n = toto+1:totn
		if (n==1)&&~flag %If this is the first iteration, it presents the user with the list of available choices
            [Names,Vals,iI] = readSETfile([SPpath filesep SPfile{n}]);
        else    %Otherwise is reuses the previous list.
            [Names,TVals,iI] = readSETfile([SPpath filesep SPfile{n}],iI);
            for m = 1:length(iI)
                Vals{m} = [Vals{m}; TVals{m}];  %Appends the information from the nth file to the list(s).
            end
            clear TVals
		end
    end
    if totn >= N    %Checks to see if another iteration is needed.
        done2 = true;
        N = totn;
    else
        toto = totn;
    end
end
