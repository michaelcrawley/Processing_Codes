function [Names,Vals] = gatherSETdata(N,varargin)
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
%
%OUTPUTS
% Names - a cell array containing the names of the data entries gathered by
% the program.
% Vals - a cell array containing the composited entries gathered by the
% program.
%
%This function can be called in two forms:
% 1) [Names,Vals] = gatherSETdata(N);  -If no starting directory is given,
% the program uses the default starting point.
% 2) [Names,Vals] = gatherSETdata(N,SD); -the program uses the user
% specificed starting directory.


if isempty(varargin)
    SD = 'My Computer';
else
    SD = varargin{1};
end

done2 = false;
toto = 0;
n = 0;
SPfile = [];
while ~done2
    
    SPpath = uigetdir(SD,'Specify Directory Containing Data');
    flist = struct2cell(dir(SPpath));
        %Eliminates all files not having .SET in name
    q = strfind(flist(1,:),'.SET');
    keep = zeros(size(q));
    for n = 1:length(keep)
        if ~isempty(q{n})
            keep(n) = 1;
        end
    end
    flist = flist(1,logical(keep));
        %Lets user choose the subset of files in that directory for gathering.
    [keep,ok] = listdlg('PromptString','Select files to process:',...
                    'SelectionMode','multiple',...
                    'ListString',flist);
    if ok==0
        error('Program Terminated Due to User Selection of Cancel')
    end
    SPf = flist(keep);
    
        %Scans through all the kept files extracting the information
    SPfile = [SPfile SPf];
    totn = toto+length(SPf);
    for n = toto+1:totn
        if n==1 %If this is the first iteration, it presents the user with the list of available choices
            [Names,Vals,kp] = readSETfile([SPpath filesep SPfile{n}]);
        else    %Otherwise is reuses the previous list.
            [TN,TVals,kp] = readSETfile([SPpath filesep SPfile{n}],kp);
            for m = 1:length(kp)
                Vals{m} = [Vals{m}; TVals{m}];  %Appends the information from the nth file to the list(s).
            end
            clear TN TVals
        end
    end
    if totn >= N    %Checks to see if another iteration is needed.
        done2 = true;
        N = totn;
    else
        toto = totn;
    end
end
