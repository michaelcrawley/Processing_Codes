function info = parse_pressure_fname(fname_w_ext,M_j,D,info_in)

%   Parsing a filename like
%       'M0.9_m(00)(00)_F0.00_PW0.0_T19.5_IA00_On00000_R00.NOS'
%   OR
%       'M0.9_m(00)(00)_F2.02_PWAUH_T16.3_IA00_On00255_R00.NOS'

%   Parsing a filename like 'M13_X01_R1p09_WF_m00_F2.10_DAUH_T12.9.lvm'.
%       Mach #, 
%       1st azimuthal mode of actuation
%       2nd azimuthal mode of actuation
%       Forcing frequency in kHz
%       Pulse width or "Auto duty cycle (AUH)"
%       Stagnation temperature in deg C
%       Increase actuator
%       Actuator on/off
%       Rotations

TC2K = 273.15; %degC to Kelvin

if nargin < 3 || isempty(D)
    D = 0.0254;
end
if nargin < 4 || isempty(info_in)
    info_in.thermo = [];
end

info = info_in;

if ~isfield(info.thermo,'gamma')
    info.thermo.gamma = 1.4; %Specific heat ratio for air
end
if ~isfield(info.thermo,'R_gas_const_air')
    info.thermo.R_gas_const_air = 287.05; %Gas constant for air, J/kg-K
end
if ~isfield(info.thermo,'P_amb')
    info.thermo.P_amb = 101325; %Standard atmospheric pressure
end
if ~isfield(info.thermo,'D')
    info.thermo.D = D; %Jet diameter
end

dot_locs = findstr(fname_w_ext,'.'); %Locate dots in filename with extension

fname = fname_w_ext(1:dot_locs(end)-1); %Get rid of extension

if ~isfield(info.thermo,'M_j')
    info.thermo.M_j = parse_local_fnc(fname,'M');
    if isempty(info.thermo.M_j) 
        if nargin > 1 && ~isempty(M_j)
            info.thermo.M_j = M_j;
        else
            info.thermo.M_j = 0;
        end
    end
end
[info.actuation.actuation_mode1,info.actuation.actuation_mode2] = parse_local_azmode(fname,'m');
info.actuation.actuation_dc = parse_local_pw(fname,'PW');
info.actuation.actuation_phase = parse_local_fnc(fname,'p');
info.thermo.T_0 = TC2K+parse_local_fnc(fname,'T');
info.actuation.increase_actuator = parse_local_fnc(fname,'IA');
info.actuation.actuator_on_off = parse_local_fnc(fname,'On');
info.actuation.rotation = parse_local_fnc(fname,'R');

%   Either the Strouhal number or the forcing frequency (in kHz) will be
%   provided
info.actuation.St_DF = parse_local_fnc(fname,'S');
info.actuation.actuation_freq = parse_local_fnc(fname,'F')*1000;

info.thermo.P_j = info.thermo.P_amb;
info.thermo.T_j = info.thermo.T_0/(1+((info.thermo.gamma-1)/2)*info.thermo.M_j^2);
info.thermo.spd_snd_j = sqrt(info.thermo.gamma*info.thermo.R_gas_const_air*info.thermo.T_j);
info.thermo.rho_j = info.thermo.P_j/(info.thermo.R_gas_const_air*info.thermo.T_j);
info.thermo.mu_j = (0.00001827*(291.15+120)/(info.thermo.T_j+120)*(info.thermo.T_j/291.15)^(3/2));
info.thermo.U_j = info.thermo.M_j*info.thermo.spd_snd_j;

if ~isfield(info.thermo,'Re')
    info.thermo.Re = info.thermo.rho_j*info.thermo.U_j*info.thermo.D/info.thermo.mu_j;
end

%   Convert from St_DF to F_f, or vice versa, depending on which one was
%   supplied
if ~isempty(info.actuation.St_DF)
    info.actuation.actuation_freq = info.actuation.St_DF*info.thermo.U_j/info.thermo.D;
elseif ~isempty(info.actuation.actuation_freq)
    info.actuation.St_DF = info.actuation.actuation_freq*info.thermo.D/info.thermo.U_j;
end

if isempty(info.actuation.St_DF) || info.actuation.St_DF == 0
    info.baseline = 1;
else
    info.baseline = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = parse_local_fnc(fname,str)

start_loc = findstr(fname,str);
if isempty(start_loc)
    val = [];
    return;
end
fname_rem = fname(start_loc+length(str):end);
ends_loc_temp = findstr(fname_rem,'_');
if isempty(ends_loc_temp)
    val = str2double(fname_rem);
else
    val = str2double(fname_rem(1:ends_loc_temp(1)-1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mode1,mode2] = parse_local_azmode(fname,str)

str = [str,'('];
mode1 = [];
mode2 = [];
start_loc = findstr(fname,str);
if isempty(start_loc)
    return;
end
fname_rem = fname(start_loc+length(str):end);
ends_loc_temp = findstr(fname_rem,')');
if isempty(ends_loc_temp)
    return;
else
    mode1 = str2double(fname_rem(1:ends_loc_temp(1)-1));
end
if fname_rem(ends_loc_temp(1)+1) ~= '('
    return;
end
fname_rem = fname_rem(ends_loc_temp(1)+2:end);
ends_loc_temp = findstr(fname_rem,')');
if isempty(ends_loc_temp)
    return;
else
    mode2 = str2double(fname_rem(1:ends_loc_temp(1)-1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = parse_local_pw(fname,str)

start_loc = findstr(fname,str);
if isempty(start_loc)
    val = -1;
    return;
end
fname_rem = fname(start_loc+length(str):end);
ends_loc_temp = findstr(fname_rem,'_');
if isempty(ends_loc_temp)
    dc_str = fname_rem;
else
    dc_str = fname_rem(1:ends_loc_temp(1)-1);
end

if strcmp(dc_str,'AUH')
    val = -1;
else
    val = str2double(dc_str);
end
