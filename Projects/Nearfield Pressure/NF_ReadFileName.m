function [phys] = NF_ReadFileName(filename,pp)
%Parses specified filename to determine jet operating conditions (TTR, exit
%velocity, acoustic Mach number, forcing frequency, etc).  Fields are added
%to the 'phys' structure.
    
    md = getMetaFromStr(filename,'NearFieldParams'); %pull metadata from filename
    
    %pull necessary data
    M = md.M.value;
    T = md.T.value;
    x = (md.x.value:(pp.sp*cosd(md.a.value)):md.x.value+(pp.nm-1)*pp.sp*cosd(md.a.value)); %axial distances of microphones
    y = md.r.value+(x-x(1))*tand(md.a.value); % y/D
    r = [pp.R sqrt(x.^2+y.^2)*0.0254];%radial distance of nearfield microphones, converted to meters
    gain = [md.mVf.value*ones(length(pp.FFCh),1); 0]; %get farfield microphone gains, add in zero for trigger.  Assumes farfield mics were at same gain, and trigger is between farfield and nearfield mics
    if isfield(md,'mV')
        gain = [gain; md.mV.value*ones(pp.nm,1)];
    elseif isfield(md,'mVu') && isfield(md,'mVd');
        gain = [gain; md.mVu.value*ones(pp.nm/2,1); md.mVd.value*ones(pp.nm/2,1)]; %assumes gain settings were split in half over the NF array
    end
    
    
    %Calc jet/ambient properties
    To = T + 273.15; %convert to kelvin
    Te = To/(1+1/5*M^2); %jet exit temperature, in kelvin
    TTR = To/(pp.AT+273.15); %jet temperature ratio
    Ue = M*sqrt(1.4*287.05*Te); %nozzle exit velocity (m/s)
    Ma = M*sqrt(Te/(pp.AT+273.15)); %acoustic jet Mach number
    a = sqrt(1.4*287.05*(pp.AT+273.15)); %ambient speed of sound (m/s)
    c = Ue/M; %speed of sound in jet (m/s)
    Uc = Ue*a/(a+c); %estimated convective velocity
    t_c = r/Uc; %convective time for structures
    i_c = round(t_c*pp.FS); %convective index for structures
    t_a = r/a; %acoustic time for actuator pulse
    i_a = round(t_a*pp.FS); %acoustic index for actuator pulse
    
    phys = struct('M',M,'To',To,'x',x,'y',y,'r',r,'TTR',TTR,'Ue',Ue,'Ma',Ma,'a',a,'c',c,'Uc',Uc,'t_c',t_c,'i_c',i_c,'t_a',t_a,'i_a',i_a,'gain',gain);
end