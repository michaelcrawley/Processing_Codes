function [pp] = ProcessCalFiles(pp)
%Reads in calibration files, processes them for each channel, and updates
%the processing parameters structure with the results.  Also determines the
%file format (ascii, binary, etc) of the data files to be processed.

%Last updated by Michael Crawley on 2011-07-04

    pp.Nch = length(pp.CalFiles); %Determines the number of channels of data by the number of calibration files
    ChN = zeros(1,pp.Nch);
    for n = 1:pp.Nch   %Extracts channel number from from file name
        tmp = regexpi(pp.CalFiles{n},'_CH');
        if isempty(tmp)
            error('Channel Number Prefix Must be CH, Ch, or ch')
        end
        tmp2 = pp.CalFiles{n}(tmp+3:tmp+8);
        tmp3 = strfind(tmp2,'_');
        ChN(n) = str2double(tmp2(1:tmp3(1)-1));
    end
    if length(unique(ChN))~=length(ChN)
        error('Calibration channel numbers are not unique')
    end
    [~,IX] = sort(ChN); clear n ChN B tmp tmp2 tmp3;
    pp.CalFiles = pp.CalFiles(IX); clear IX; %Reorders files in ascending channel number to ensure proper correlation to data columns
    
    pp.ChRef = zeros(1,pp.Nch);
    for n = 1:pp.Nch   %Calculates spectral reference value for each channel from the calibration files
        if ~exist('isbinary','var')
            isbinary = true;
            fid = fopen([pp.CalPath '\' pp.CalFiles{n}],'r');
            tmp = fread(fid,500,'*char'); fclose(fid);
            if strcmp('\t',pp.delimI)
                if length(strfind(tmp',char(9))) > 5
                    isbinary = false;
                end
            else
                if length(strfind(tmp',pp.delimI)) > 5
                    isbinary = false;
                end
            end
            if isbinary
                if length(strfind(tmp',char(10))) > 5	%checks for line feed characters if delimiter wasn't found
                    isbinary = false;
                end
            end
        end
        if isbinary
            fid = fopen([pp.CalPath '\' pp.CalFiles{n}],'r');
            tmp = fread(fid,'float32'); fclose(fid);

            Ltmp = length(tmp);
            if round(Ltmp/pp.BS/pp.Nch)==Ltmp/pp.BS/pp.Nch	%If there could be more than one channel, unpack as if there are "Nch" channels
                if pp.TS_cal	%If timestamp is present in files
                    tmp2 = reshape(tmp,pp.BS,pp.Nch+1,[]);	%parses data into channels along second dimension and blocks along third dimension
                    tmp2 = permute(tmp2,[1 3 2]);	%reorders data into channels along third and blocks along second dimension
                    tmp2 = reshape(tmp2,[],pp.Nch+1);	%reshapes data into 2-D matrix with all blocks for a channel in one column
                    tmp2 = tmp2(:,2:end);	%removes timestamp
                else
                    tmp2 = reshape(tmp,pp.BS,pp.Nch,[]);
                    tmp2 = permute(tmp2,[1 3 2]);
                    tmp2 = reshape(tmp2,[],pp.Nch);
                end
                mtmp = mean(tmp2.^2);	%If it was actually only one channel, all imagined channels will have same mean square. In which case the sum below will equal the number of channles. 
                if sum(round(mtmp/mtmp(n))) ~= pp.Nch
                    tmp = tmp2(:,n);
                    disp(['     Cal File ' num2str(n) ' has ' num2str(pp.Nch) ' channels in it']);
                end
                clear tmp2 mtmp;
            else
                disp(['     Cal File ' num2str(n) ' has 1 channel in it']);
            end
        else
            tmp = dlmread([pp.CalPath '\' pp.CalFiles{n}],pp.delimI);
            if size(tmp,2)~=1   %If there is only one channel of data, assume it is for the nth channel    
                %Else, assume that the nth column contains the calibration for the nth channel
                tmp = tmp(:,n);
            end
        end
        S = length(tmp);
        if round(S/pp.BS)~=S/pp.BS
            error(['Data in ' pp.CalFiles{n} ' is not an integer number of blocks'])
        end
        tmp = reshape(tmp,pp.BS,[]);    %Cut the data into blocks
        tmp = std(tmp,0,1);     %Remove mean and compute RMS of blocks
        pp.ChRef(n) = mean(tmp.^2);    %Compute mean square of prms of blocks (i.e. prms^2)
    end
    pp.ChRef = spdiags(1./pp.ChRef(:),0,pp.Nch,pp.Nch); %prepares the calibration values for scaling the data.

        %Inserts extra \ into directory path so it will print properly
    tmp = fliplr(strfind(pp.CalPath,'\')); pp.CalPathP = pp.CalPath;
    for n = 1:length(tmp)
        pp.CalPathP = [pp.CalPathP(1:tmp(n)) pp.CalPathP(tmp(n):end)];
    end
    disp('  Calibration Data Processed...Beginning Data Processing')
end