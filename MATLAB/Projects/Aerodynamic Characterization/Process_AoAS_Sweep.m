function Process_AoAS_Sweep(src,params)

%%%%%%%%%% Necessary fields for params file
%   P_inf           ambient pressure (pascals)
%   Temp            ambient temperature (kelvin)
%   Velocity        Wind tunnel velocity (m/s)
%   loadcell        type of load cell used, and orientation
%   Area            Reference Surface Area
%   Chord           Reference Chord length
%   span            Reference Span length
%   Xmrp            X MRP from sensor to body, Positive Downstream
%   Ymrp            Y MRP from Sensor to body, Positive Out Right wing
%   Zmrp            Z MRP from Sensor to body, Positive Downward
%   AoAB            Angle of attack offset for balance
%   AoAM            Angle of attack offset for model
%   AoSB            Sideslip of balance from flow
%   AoSM            Sideslip of Model from flow
%   AoRB            Roll offset of balance
%   AoRM            Roll offset of model

    %Constants
    Rot1 = @(angle) [cos(angle) sin(angle) 0;  -sin(angle) cos(angle) 0; 0 0 1];
    Rot2 = @(angle) [cos(angle) 0 sin(angle) ; 0 1 0 ; -sin(angle) 0 cos(angle)];
    Rot3 = @(angle) [1 0 0 ; 0 cos(angle) sin(angle) ; 0 -sin(angle) cos(angle)];

    %%% Get Calibration Data
    [cal,f_indx,T_indx] = Calibration(params.loadcell);
    
    %%% Find Data and Tare files
    tare_files = getfiles('Tare*_Raw.wtd',[src,filesep,'Tare']);
    num_tares = length(tare_files);
    data_files = getfiles('*_Raw.wtd',src);
    num_data = length(data_files);
    if num_tares ~= num_data
        error('Mismatch between number of tares and number of data files!')
    end
    
    %%% Process constant parameters
    rho = params.P_inf/(params.Temp*287.05);          %Calculated Gas Constant 
    Q = 0.5 * rho * params.Velocity^2;                           %Calculated Dynamic Pressure
    AoSB = params.AoSB*pi/180;                         %convert Sideslip of Balance from flow to radians
    AoSM = (params.AoSB + params.AoSM)*pi/180;      %Sideslip of Model from flow
    RollB = params.AoRB*pi/180;                       %Roll of Balance from wings level flow
    RollM = (params.AoRB + params.AoRM)*pi/180;    %Roll of Model from wings level flow

    delAoS = AoSM-AoSB;    %Angle difference from X axis
    delRoll = RollM-RollB; %Angle difference from z axis
    
    %%% Process Wind-Off Tares
    tare(1:num_tares) = struct();
    for n = 1:num_tares
        %Get Angle of Attack & Sideslip
        angles = regexp(tare_files{n},'_AoA(?<AoA>[0-9.-]+)+_AoS(?<AoS>[0-9.-]+)','names');
        tare(n).AoA = str2double(angles.AoA);
        tare(n).AoS = str2double(angles.AoS);
        
        %Read voltages & convert to mean load
        fid = fopen([src,filesep,'Tare',filesep,tare_files{n}]);
        tmp = textscan(fid,'%f','headerlines',1);
        fclose(fid);
        voltage = reshape(tmp{1},6,[]);
        tare(n).means = mean(cal*voltage,2)';        
    end
    %sort by increasing angle of attack
    [~,indx] = sort([tare.AoA]); 
    tare = tare(indx);
    
    %%% Process Loads Data
    data(1:num_data) = struct();
    for n = 1:num_data
        %Get Angle of Attack & Sideslip
        angles = regexp(data_files{n},'_AoA(?<AoA>[0-9.-]+)+_AoS(?<AoS>[0-9.-]+)','names');
        data(n).AoA = str2double(angles.AoA);
        data(n).AoS = str2double(angles.AoS);
        
        %Read raw voltages & convert to force/torque
        fid = fopen([src,filesep,data_files{n}]);
        tmp = textscan(fid,'%f','headerlines',1);
        fclose(fid);
        voltage = reshape(tmp{1},6,[]);
        loads = (cal*voltage)';
        data(n).raw.mean = mean(loads,1);
        data(n).raw.std = std(loads,[],1);
        
        %Interp between AoA's in tares
        int_tare = interp1([tare.AoA],reshape([tare.means],6,num_tares)',data(n).AoA,'linear','extrap');
        data(n).raw.corrected = loads - repmat(int_tare,size(loads,1),1);
        data(n).raw.corrected_mean = mean(data(n).raw.corrected,1);
        
        %Calculate Rotation Matrices
        AoAB = (data(n).AoA+params.AoAB)*pi/180;                %AoA offset of balance from inclinometer
        AoAM = (data(n).AoA+params.AoAB+params.AoAM)*pi/180;   %AoA offset of model from balance
        delAoA = AoAM-AoAB;    %Angle difference from y axis
        aircraft_body_rot = Rot3(delRoll)*Rot2(delAoA)*Rot1(delAoS);
        stability_rot = Rot1(AoSM) * Rot2(AoAM);
        
        %Convert axes on f/T
        forces = f_indx(data(n).raw.corrected)';
        moments = T_indx(data(n).raw.corrected)';
        forces = stability_rot*aircraft_body_rot*forces;
        moments = stability_rot*aircraft_body_rot*moments;
        MRP = stability_rot*aircraft_body_rot*[params.Xmrp ; params.Ymrp ; params.Zmrp];
        
        D = forces(1,:)';
        S = forces(2,:)';
        L = forces(3,:)';
        RMT = moments(1,:)';
        PMT = moments(2,:)';
        YMT = moments(3,:)';
        XP = MRP(1);
        YP = MRP(2);
        ZP = MRP(3);
        
        %Calculate Aerodynamic coefficients
        RM = RMT+L*YP-S*ZP;
        PM = PMT-L*XP-D*ZP;
        YM = YMT-D*YP-S*XP;

        CL = L/Q/params.Area;
        CD = D/Q/params.Area;
        CS = S/Q/params.Area;

        CRM = RM/Q/params.Area/params.span;
        CPM = PM/Q/params.Area/params.Chord;
        CYM = YM/Q/params.Area/params.span;
        
        data(n).aero.CL = mean(CL);
        data(n).aero.CD = mean(CD);
        data(n).aero.CS = mean(CS);
        data(n).aero.CRM = mean(CRM);
        data(n).aero.CPM = mean(CPM);
        data(n).aero.CYM = mean(CYM);
        
        data(n).aero.CLstd = std(CL);
        data(n).aero.CDstd = std(CD);
        data(n).aero.CSstd = std(CS);
        data(n).aero.CRMstd = std(CRM);
        data(n).aero.CPMstd = std(CPM);
        data(n).aero.CYMstd = std(CYM);
    end
    %sort by increasing angle of attack
    [~,indx] = sort([data.AoA]); 
    data = data(indx); %#ok<NASGU>

    processing_file = [mfilename '.m']; %#ok<NASGU>
    processing_date = date; %#ok<NASGU>
    save([src filesep 'Results.mat'],'params','src','data','tare','processing_file','processing_date');
end

function [cal,f_indx,T_indx] = Calibration(loadcell)
    switch lower(loadcell)
        case 'jr3_load_cell' 
            cal = importdata('JR3_CalibrationMatrix.txt');
            f_indx = @(x) x(:,[3,2,1]).*repmat([-1 1 -1],size(x,1),1);
            T_indx = @(x) x(:,[6,5,4]).*repmat([-1 1 1],size(x,1),1);
        case 'jr3_fz_up'
            cal = importdata('JR3_CalibrationMatrix.txt');
            f_indx = @(x) x(:,[1,2,3]).*repmat([1 -1 1],size(x,1),1);
            T_indx = @(x) x(:,[4,5,6]).*repmat([1 -1 1],size(x,1),1);
        case 'ati_n25_fz_up'
            cal = importdata('ATI_N25_FT14574.txt');
            f_indx = @(x) x(:,[1,2,3]);
            T_indx = @(x) x(:,[4,5,6]);
        case 'ati_n25_fz_down'
            cal = importdata('ATI_N25_FT14574.txt');
            f_indx = @(x) x(:,[1,2,3]).*repmat([1 -1 -1],size(x,1),1);
            T_indx = @(x) x(:,[4,5,6]).*repmat([1 -1 -1],size(x,1),1);
        otherwise
            error('Incorrect load cell definition');
    end
end