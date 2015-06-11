function [Names,Vals,keep] = readSETfile(fname,keep)
%This function reads the contents of a SET file. In order to operate
%properly, the SET file must be formatted as a tab delimited file in which
%the first column contains the data label for a given row.
%
%INPUT
% fname - the full path to the SET file to be read.
% keep - a vector containing the line numbers of the data items which should be outputted.
%
%OUTPUT
% Names - a cell array containing the names of the data entries gathered by
% the program.
% Vals - a cell array containing the composited entries gathered by the
% program.
% keep - a vector containing the rows in the set file which were outputted. 
%
%There are two valid calls to this function:
% 1) [Names,Vals,keep] = readSETfile(fname); -Since there is not vector
% telling the program which entries are to be outputted, the program
% presents the user with the list of available choices.
% 2) [Names,Vals,keep] = readSETfile(fname,keep); -The program uses the
% entries in the vector "keep" to determine which entries in the SET file
% are to be outputted.

    %This section parses the file name for useful information.
[path,name] = fileparts(fname);
q = strfind(name,'_');
if ~isempty(q)
    fp{1} = name(1:q(1)-1);
    if length(q) > 1
        for n = 2:length(q)
            fp{n} = name(q(n-1)+1:q(n)-1);
        end
    else
        n = 1;
    end
    fp{n+1} = name(q(end)+1:end);
    
    k = 0;
    for n = 1:length(fp)
        k = k+1;
        if strmatch(fp{n}(1),'M')   %Mach number is already in the SET file
            k = k-1;
        elseif strmatch(fp{n}(1),'F')   %Forcing Frequency is already in the SET file
            k = k-1;
		elseif strmatch(fp{n}(1),'S')   %Forcing Strouhal number is already in the SET file
            k = k-1;
        elseif strmatch(fp{n}(1),'T')   %Stagnation temperature is already in the SET file
            k = k-1;
        elseif strmatch(fp{n}(1:2),'Ba')  %Ignores Baseline tag
            k = k-1;
        elseif strmatch(fp{n}(1),'m')
            Names{k} = 'Azimuthal Mode';
            Vals{k} = fp{n}(2:end);
        elseif strmatch(fp{n}(1:2),'PW')
            Names{k} = 'Pulse Width';
            Vals{k} = fp{n}(3:end);
        elseif strmatch(fp{n}(1:2),'IA')
            Names{k} = 'Increment Actuators';
            Vals{k} = str2num(fp{n}(3:end));
        elseif strmatch(fp{n}(1:2),'On')
            Names{k} = 'Actuators On';
            Vals{k} = str2num(fp{n}(3:end));
        elseif strmatch(fp{n}(1),'R')
            Names{k} = 'Pulses per Drift';
            Vals{k} = str2num(fp{n}(2:end));
        else
            Names{k} = fp{n};
            Vals{k} = fp{n};
        end
    end
end

    %This section goes through the file one line at a time, parsing the
    %information into Names and Values.
fid = fopen(fname,'r');
done = false;   n = k;  good = ones(1,k);
q = fgetl(fid);
while ~done
    sep = regexp(q,'\t');
    n = n+1;
    if ~isempty(sep)    %Only parses entries that are tab delimited
        good(n) = 1;
        Names{n} = q(1:sep(1)-1);   %The name is always the first column
        if length(sep)==1   %If there is only one data value in the row
            qq = q(sep+1:end);
            w = unique(qq);
            w = w(isletter(w));
            if isempty(w) || ((length(w)==1) && strcmpi(w,'e'))   %If the values are numeric, convert to numbers
                errSym = strfind(qq,'+/-'); %Check for the plus-minus symbol present in uncertainty estimates of slopes, etc.
                if isempty(errSym)
                    Vals{n} = str2num(qq);
                else
                    Vals{n} = [str2num(qq(1:errSym-2)) str2num(qq(errSym+4:end))];
                end
            else
                Vals{n} = qq;
            end
        else    %If there are multiple data values in the row, iterate through them all
            clear d
            for m = 1:length(sep)
                if m < length(sep)
                    qq = q(sep(m)+1:sep(m+1)-1);
                else
                    qq = q(sep(m)+1:end);
                end
                w = unique(qq);
                w = w(isletter(w));
                if isempty(w) || ((length(w)==1) && strcmpi(w,'e'))   %If the values are numeric, convert to numbers
                    errSym = strfind(qq,'+/-');
                    if isempty(errSym)
                        d{m} = str2num(qq);
                    else
                        d{m} = [str2num(qq(1:errSym-2)) str2num(qq(errSym+4:end))];
                    end
                else
                    d{m} = qq;
                end
            end
            Vals{n} = d;
        end
    else
        good(n) = 0;
    end
    q = fgetl(fid);
    if q == -1  %reads until end of file is encountered
        done = true;
    end
end
fclose(fid);
good = logical(good);

    %If keep is not already specified by user input, present the user with
    %the available options.
if ~exist('keep','var')
    [keep,ok] = listdlg('PromptString','Data Columns:',...
                    'SelectionMode','multiple',...
                    'ListString',Names(good));
    if ok==0
        error('Program Terminated Due to User Selection of Cancel')
    end
    ng = find(~good);
    for n = 1:length(keep)
        keep(n) = keep(n) +sum(keep(n) > ng);
    end
end
Names = Names(keep);    %Keep only the user selected options
Vals = Vals(keep);

    %Scan through each kept entry to see if the information should be
    %converted from a cell array to a matrix.
for n = 1:length(Vals)  
    if iscell(Vals{n})
        if isnumeric(Vals{n}{1})
            Vals{n} = cell2mat(Vals{n});
        end
    end
end
