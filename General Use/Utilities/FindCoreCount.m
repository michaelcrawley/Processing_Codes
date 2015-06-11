function corenum = FindCoreCount
    [status, tmpstr] = dos('wmic cpu get NumberOfLogicalProcessors'); %verified on Windows 7 x64
    
    if status == 0
       corenum = str2double(tmpstr(26:end)); %parses resulting string
    else
        error('Unable to find number of cores');
    end
end