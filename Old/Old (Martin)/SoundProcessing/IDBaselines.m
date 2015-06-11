function flistn = IDBaselines(fdir,flist)
%This function scans through the list of .NOS files looking for entries
%with zero forcing frequency. These zero frequency entries are stripped of
%their extraneous information and "_Baseline" is added immediately following
%the Mach number. Note that this will not remove a suffix or any other
%element that doesn't conform to the file naming standards. Typical result
%looks like: 'Mxx_Baseline_Txxxx.NOS'.

flistn = flist;
[fp,fn,EXT] = fileparts(flist{1}); clear fp fn;
fix = 0;
for n = 1:length(flist)
    Fnm = flist{n};
    FLoc = strfind(Fnm,'_F');	%Looks for forcing frequency present in file name
	SLoc = strfind(Fnm,'_S');	%Looks for forcing Strouhal number present in file name
    if ~isempty(FLoc) || ~isempty(SLoc)
        if isempty(FLoc)
			Fend = min(strfind(Fnm(SLoc+1:end),'_'))+SLoc;
			FH = str2num(Fnm(SLoc+2:Fend-1)); clear FLoc SLoc Fend;
		else
			Fend = min(strfind(Fnm(FLoc+1:end),'_'))+FLoc;
			FH = str2num(Fnm(FLoc+2:Fend-1))*1000; clear FLoc SLoc Fend;
		end
        if FH==0	%If forcing frequency is zero, strip forcing info from name
			N = Fnm; tmp = length(N)-length(EXT);
			N = [N(1:tmp) '_' N(tmp+1:end)];	%Extra underscore temporarily added to file name before extension for ease of stripping
            uScore = strfind(N,'_');
			for m = length(uScore)-1:-1:1
				if strcmp(N(uScore(m)+1),'m')	%Strip azimuthal mode
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];
				elseif strcmp(N(uScore(m)+1),'F')	%Strip forcing frequency
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];
				elseif strcmp(N(uScore(m)+1),'S')	%Strip forcing Strouhal number
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];
				elseif strcmp(N(uScore(m)+1:uScore(m)+2),'PW')	%Strip pulse width
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];	
				elseif strcmp(N(uScore(m)+1:uScore(m)+2),'IA')	%Strip increment actuators
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];
				elseif strcmp(N(uScore(m)+1:uScore(m)+2),'On')	%Strip actuators on
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];
				elseif strcmp(N(uScore(m)+1),'R')	%Strip rotate actuators
					N = [N(1:uScore(m)-1) N(uScore(m+1):end)];
				end
			end
			tmp = length(N)-length(EXT);
			N = [N(1:uScore(1)-1) '_Baseline' N(uScore(1):tmp-1) EXT];	%Add "_Baseline" to filename and remove temporary underscore
				
			TLoc = strfind(N,'_T');   
            Tend = min(strfind(N(TLoc+1:end),'_'))+TLoc;
            if isempty(Tend)
                Tend = length(N)-length(EXT)+1;
            end
            Ttxt = N(TLoc+2:Tend-1);
            Ttxt(strfind(Ttxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
            Tmp = str2num(Ttxt);
            
            mtch = sum(strcmp(N,flistn)); %If the resulting file has the same filename as another, add another decimal place to the temperature which will make it unique.
            if mtch ~= 0
                fix = fix +1;
                if round(Tmp)==Tmp
                    TF = ['.' num2str(fix)];
                else
                    TF = num2str(fix);
                end
                N = [N(1:TLoc+1) Ttxt TF N(Tend:end)];
            end
            flistn{n} = N;
            movefile(fullfile(fdir,Fnm),fullfile(fdir,N));
        end
    end
end

flistn = sort(flistn);	%Sorts list in ASCII alphabetical order
