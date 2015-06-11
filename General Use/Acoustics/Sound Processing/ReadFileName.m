function [data] = ReadFileName(filename,pp)
%Parses specified filename to determine jet operating conditions (TTR, exit
%velocity, acoustic Mach number, forcing frequency, etc).  Fields are added
%to the 'data' structure.

%Last updated by Michael Crawley on 2012-02-09

    disp(['    Processing File: ' filename]);
        %Locates temperature information in file name
    TLoc = strfind(filename,'_T');   
    Tend = min(strfind(filename(TLoc+1:end),'_'))+TLoc;
    if isempty(Tend)
        Tend = length(filename)-length(['.' pp.EXT])+1;
    end
    Ttxt = filename(TLoc+2:Tend-1); clear TLoc Tend;
    Ttxt(strfind(Ttxt,'p'))='.';    %if p was used instead of . in floating point number, it is replaced
    data.To = floor(str2num(Ttxt)*10)/10+273.15; clear Ttxt; %#ok<ST2NM> %Converts to Kelvin, truncates to first decimal place
	data.Temp = data.To/(1+1/5*pp.M^2);  %Calculates isentropically expanded jet exit temperature - K
	
		%Locates forcing frequency information in file name if it exists
	FLoc = strfind(filename,'_F');  
	SLoc = strfind(filename,'_S');
    if ~isempty(FLoc)	%Looks for Forcing frequency (kHz) in filename
        data.SHarmonics=true;    %This boolean will be used to determine if forcing frequency information needs to be written in output file
        Fend = min(strfind(filename(FLoc+1:end),'_'))+FLoc;
        data.FH = str2num(filename(FLoc+2:Fend-1))*1000; %#ok<ST2NM>
		data.StDF = data.FH*pp.D/(pp.M*sqrt(7/5*287.05*data.Temp));
	elseif ~isempty(SLoc)	%Looks for Forcing strouhal number in filename
        data.SHarmonics=true;    %This boolean will be used to determine if forcing frequency information needs to be written in output file
        Fend = min(strfind(filename(SLoc+1:end),'_'))+SLoc;
        data.StDF = str2num(filename(SLoc+2:Fend-1)); %#ok<ST2NM>
		data.FH = data.StDF*pp.M*sqrt(7/5*287.05*data.Temp)/pp.D;
	else
        data.SHarmonics=false;   
        data.FH = 0;
    end
    clear FLoc SLoc Fend;
    
    data.TTR = data.To/(pp.AT+273.15);
    data.Ue = pp.M*sqrt(1.4*287.05*data.Temp);
    data.Mj = pp.M*sqrt(data.Temp/(pp.AT+273.15));
    data.rhoe = 1.01e5/(287.05*data.Temp);
    data.mue = (1.827e-5*(0.555*291.15+66.667)/(0.555*data.Temp+66.667).*(data.Temp/291.15).^(3/2));    
    data.Std = pp.fx*pp.D/(pp.M*sqrt(7/5*287.05*data.Temp));    %Creates Strouhal number axis assuming Mach number is constant
    data.ReN = data.rhoe*data.Ue*pp.D /data.mue;  %Calculates jet exit Reynolds number using Sutherlands formula for viscosity and ideal gas law for density
end