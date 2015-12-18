function OASPL = OASPL_calc()

dir_name = uigetdir('My Computer','Specify Directory Containing Data'); 
flist = struct2cell(dir(dir_name));
q = strfind(flist(1,:),'.fftNOS'); %extracts list of data files - ignores all other files
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

for n = 1:length(keep)
    Fnm = flist{keep(n)};
    rRaw = dlmread([dir_name filesep flist{keep(n)}],'\t',1,0); %reads data file
    Std = rRaw(:,2);
    S = size(rRaw);
    NCh = S(2)-2;
    Bg = min(find(Std > 0.04)); %Excludes spectral components below 0.04 Strouhal number
    Ed = min(find(Std > 4));    %Excludes spectral components above 4 Strouhal number
    for nn = 1:NCh
        OASPL(n,nn) = 10*log10(sum(10.^(rRaw(Bg:Ed,nn+2)/10)));
    end
end