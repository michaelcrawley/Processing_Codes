function flistn = IDBaselines(fdir,flist)
%This function scans through the list of .NOS files looking for entries
%with zero forcing frequency. These zero frequency entries are renamed
%using the "Baseline" format: 'Mxx_Baseline_Txxxx.NOS'.

flistn = flist;
fix = 0;
for n = 1:length(flist)
    Fnm = flist{n};
    FLoc = strfind(Fnm,'_F');   
    if ~isempty(FLoc)
        Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
        FH = str2num(Fnm(FLoc+2:Fend-1))*1000; clear FLoc Fend;
        if FH==0
            TLoc = strfind(Fnm,'_T');   
            Tend = min(strfind(Fnm(TLoc+1:end),'_'))+TLoc;
            if isempty(Tend)
                Tend = length(Fnm)-length('.NOS')+1;
            end
            Ttxt = Fnm(TLoc+2:Tend-1);
            Ttxt(strfind(Ttxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
            
            U = min(strfind(Fnm,'_'));
            
            N = [Fnm(1:U) 'Baseline_T' Ttxt '.NOS'];
            mtch = sum(strcmp(N,flistn))-1; %If the resulting file has the same filename as another, add another decimal place to the temperature which will make it unique.
            if mtch ~= 0
                fix = fix +1;
                N = [Fnm(1:U) 'Baseline_T' Ttxt num2str(fix) '.NOS'];
            end
            flistn{n} = N;
            movefile(fullfile(fdir,Fnm),fullfile(fdir,N));
        end
    end
end

flistn = struct2cell(dir(fdir));
q = strfind(flistn(1,:),'.NOS'); %extracts list of data files - ignores all other files
keep = ~cellfun('isempty',q);
flistn = flistn(1,logical(keep));