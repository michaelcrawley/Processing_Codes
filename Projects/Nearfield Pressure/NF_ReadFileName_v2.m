function [phys] = NF_ReadFileName_v2(filename,pp)
%Parses specified filename to determine jet operating conditions (TTR, exit
%velocity, acoustic Mach number, forcing frequency, etc).  Fields are added
%to the 'phys' structure. This version will account for the different
%time-of-arrivals from each individual actuator, as well as for different
%forcing modes
    
    md = getMetaFromStr(filename,'NearFieldParams_v2'); %pull metadata from filename
    
    %Standard experimental configuration parameters
    D = 0.0254;
    N = 8;
    act_loc = 0:pi/4:7*pi/4; %actuator locations
    act_y = (D/2)*sin(act_loc);
    act_z = (D/2)*cos(act_loc);
    NF_phi = 0;
    FF_phi = 15; %in degrees
    
    %pull necessary data
    M = md.M.value;
    T = md.T.value;
    Stdf = md.S.value;
    x = (md.x.value:(pp.sp*cosd(md.a.value)):md.x.value+(pp.nm-1)*pp.sp*cosd(md.a.value)); %axial distances of microphones
    y = md.r.value+(x-x(1))*tand(md.a.value); % y/D
    r = [pp.R sqrt(x.^2+y.^2)*D];%radial distance of nearfield microphones, converted to meters
    gain = [md.mVf.value*ones(length(pp.FFCh),1); 0]; %get farfield microphone gains, add in zero for trigger.  Assumes farfield mics were at same gain, and trigger is between farfield and nearfield mics
    if isfield(md,'mV')
        gain = [gain; md.mV.value*ones(pp.nm,1)];
    elseif isfield(md,'mVu') && isfield(md,'mVd');
        gain = [gain; md.mVu.value*ones(pp.nm/2,1); md.mVd.value*ones(pp.nm/2,1)]; %assumes gain settings were split in half over the NF array
    end
    
    %Parse forcing mode, get forcing phases
    mode1 = str2double(md.m.value(2:3));
    mode2 = str2double(md.m.value(6:7));
    IN = 0:N-1;
    Sp1 = rem(2*pi*mode1*IN/N,2*pi);
    Sp2 = rem(2*pi*mode2*IN/N,2*pi);
    AN = 1:N;
    Sp1(Sp1<0) = Sp1(Sp1<0)+2*pi;
    Sp2(Sp2<0) = Sp2(Sp2<0)+2*pi;
    Q = abs(Sp1-Sp2) > pi;
    P = (Sp1+Sp2)/2;
    A = abs(Sp1-Sp2) <=  11*pi/16;
    P(Q == 1) = P(Q == 1)-pi;
    A(Q == 1) = 2*pi-abs(Sp1-Sp2) <=  11*pi/16;
        
    
    
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
    
    phys = struct('M',M,'To',To,'x',x,'y',y,'r',r,'TTR',TTR,'Ue',Ue,'Ma',Ma,'a',a,'c',c,'Uc',Uc,'t_c',t_c,'i_c',i_c,'t_a',t_a,'i_a',i_a,'gain',gain,'Stdf',Stdf);
end