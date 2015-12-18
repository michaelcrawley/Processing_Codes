T = '015_';

dir_name = uigetdir('My Computer','Specify Directory Containing Data'); 
flist = struct2cell(dir(dir_name));
q = strfind(flist(1,:),'.DATA'); %extracts list of data files - ignores all other files
keep = zeros(size(q));
for n = 1:length(keep)
    if ~isempty(q{n})
        keep(n) = 1;
    end
end
flist = flist(1,logical(keep));
[keep,ok] = listdlg('PromptString','Select files to process:',...
                'SelectionMode','multiple',...
                'ListString',flist);    %Asks the user to select the subset of data files they wish to process
if ok==0
    error('PT_UserCancel','Program Terminated Due to User Selection of Cancel')
end
flist = flist(keep);

for n = 1:length(keep)
    movefile([dir_name filesep flist{n}],[dir_name filesep flist{n}(1:7) T flist{n}(11:end)]);
end