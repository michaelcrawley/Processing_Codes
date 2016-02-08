function aero = Process_AoAS_Sweep(src,params)

%%%%%%%%%% Necessary fields for params struct
%   Velocity        Wind tunnel velocity (m/s)
%   loadcell        type of load cell used, and orientation
%   Area            Reference Surface Area
%   Chord           Reference Chord length
%   span            Reference Span length
%%%%%%%%%% Optional fields for the params struct
%   Xmrp            X MRP from sensor to body, Positive Downstream
%                   (default=0)
%   Ymrp            Y MRP from Sensor to body, Positive Out Right wing
%                   (default=0)
%   Zmrp            Z MRP from Sensor to body, Positive Downward
%                   (default=0)
%   AoAB            Angle of attack offset for balance
%                   (default=0)
%   AoAM            Angle of attack offset for model
%                   (default=0)
%   AoSB            Sideslip of balance from flow
%                   (default=0)
%   AoSM            Sideslip of Model from flow
%                   (default=0)
%   AoRB            Roll offset of balance
%                   (default=0)
%   AoRM            Roll offset of model
%                   (default=0)
%%%%%%%%%%%%
%last edited by mcrawley on 2016-01-22

    %Set defaults, if they do not exist
    if ~exist('params','var'), params = struct; end
    if ~isfield(params,'Velocity'), params.Velocity = input('Error: Missing Parameters! Please input freestream Velocity [m/s]: '); end
    if ~isfield(params,'loadcell')
        cell_choice = input('Error: Missing Parameters! Please choose loadcell/orientation: \n [1]: JR3 - Fz up \n [2]: JR3 - Fz down \n [3]: ATI Nano25 - Fz up \n [4]: ATI Nano25 - Fz down \n [?]:');
        switch cell_choice
            case 1
                params.loadcell = 'jr3_fz_up';
            case 2
                params.loadcell = 'jr3_fz_down';
            case 3
                params.loadcell = 'ati_n25_fz_up';
            case 4
                params.loadcell = 'ati_n25_fz_down';
        end
    end
    if ~isfield(params,'Area'), params.Area = input('Error: Missing Parameters! Please input airfoil area [m^2]: '); end
    if ~isfield(params,'Chord'), params.Chord = input('Error: Missing Parameters! Please input airfoil chord [m]: '); end
    if ~isfield(params,'span'), params.span = input('Error: Missing Parameters! Please input airfoil span [m]: '); end
    if ~isfield(params,'Xmrp'), params.Xmrp = 0; end
    if ~isfield(params,'Ymrp'), params.Ymrp = 0; end
    if ~isfield(params,'Zmrp'), params.Zmrp = 0; end
    if ~isfield(params,'AoAB'), params.AoAB = 0; end
    if ~isfield(params,'AoAM'), params.AoAM = 0; end
    if ~isfield(params,'AoSB'), params.AoSB = 0; end
    if ~isfield(params,'AoSM'), params.AoSM = 0; end
    if ~isfield(params,'AoRB'), params.AoRB = 0; end
    if ~isfield(params,'AoRM'), params.AoRM = 0; end

    %Constants
    Rot1 = @(angle) [cos(angle) sin(angle) 0;  -sin(angle) cos(angle) 0; 0 0 1];
    Rot2 = @(angle) [cos(angle) 0 sin(angle) ; 0 1 0 ; -sin(angle) 0 cos(angle)];
    Rot3 = @(angle) [1 0 0 ; 0 cos(angle) sin(angle) ; 0 -sin(angle) cos(angle)];

    %%% Get Calibration Data
    [cal,f_indx,T_indx] = Calibration(params.loadcell);

    %%% Find Data and Tare files (and figure out versioning)
    tare_files = getfiles('*.tare',src);
    num_tares = length(tare_files);
    data_files = getfiles('*.raw',src);
    num_data = length(data_files);
    version = 9; %we are assuming that the data files are written in binary
    tare_dir = src;
    if num_data == 0 %meaning we are using the old text-file based aero VI
        tare_files = getfiles('Tare*_Raw.wtd',[src,filesep,'Tare']);
        num_tares = length(tare_files);
        data_files = getfiles('*_Raw.wtd',src);
        num_data = length(data_files);
        tare_dir = [src,filesep,'Tare'];
        version = 8;
    end
    if num_tares ~= num_data
        error('Mismatch between number of tares and number of data files!')
    end

    %%% Get tunnel conditions
    tmp = regexp(data_files{1},'_','split');
    mainfname = tmp{1};
    fid = fopen([src filesep mainfname '.wtd'],'r');
    fgetl(fid);
    tmp = fscanf(fid,'%f');
    params.P_inf = tmp(1);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    tmp = fscanf(fid,'%f');
    params.Temp = tmp(3);
    fclose(fid);

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
        fid = fopen([tare_dir,filesep,tare_files{n}]);
        if version == 8
            tmp = textscan(fid,'%f','headerlines',1);
            voltage = reshape(tmp{1},6,[]);
        else
            tmp = fread(fid,'float32');
            voltage = reshape(tmp,[],6)';
        end
        fclose(fid);
        tare(n).means = mean((cal*voltage)',1);
    end
    %sort by increasing angle of attack
    [~,indx] = sort([tare.AoA]);
    tare = tare(indx);

    %%% Process Loads Data
    aero(1:num_data) = struct();
    for n = 1:num_data
        %Get Angle of Attack & Sideslip
        angles = regexp(data_files{n},'_AoA(?<AoA>[0-9.-]+)+_AoS(?<AoS>[0-9.-]+)','names');
        aero(n).AoA = str2double(angles.AoA);
        aero(n).AoS = str2double(angles.AoS);

        %Read raw voltages & convert to force/torque
        fid = fopen([src,filesep,data_files{n}]);
        if version == 8
            tmp = textscan(fid,'%f','headerlines',1);
            voltage = reshape(tmp{1},6,[]);
        else
            tmp = fread(fid,'float32');
            voltage = reshape(tmp,[],6)';
        end
        fclose(fid);
        loads = (cal*voltage)';
        raw.mean = mean(loads,1);
        raw.std = std(loads,[],1);

        %Interp between AoA's in tares
        int_tare = interp1([tare.AoA],reshape([tare.means],6,num_tares)',aero(n).AoA,'linear','extrap');
        raw.corrected = loads - repmat(int_tare,size(loads,1),1);
        raw.corrected_mean = mean(raw.corrected,1);

        %Calculate Rotation Matrices
        AoAB = (aero(n).AoA+params.AoAB)*pi/180;                %AoA offset of balance from inclinometer
        AoAM = (aero(n).AoA+params.AoAB+params.AoAM)*pi/180;   %AoA offset of model from balance
        delAoA = AoAM-AoAB;    %Angle difference from y axis
        aircraft_body_rot = Rot3(delRoll)*Rot2(delAoA)*Rot1(delAoS);
        stability_rot = Rot1(AoSM) * Rot2(AoAM);

        %Convert axes on f/T
        forces = f_indx(raw.corrected)';
        moments = T_indx(raw.corrected)';
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

        aero(n).CL = mean(CL);
        aero(n).CD = mean(CD);
        aero(n).CS = mean(CS);
        aero(n).CRM = mean(CRM);
        aero(n).CPM = mean(CPM);
        aero(n).CYM = mean(CYM);

        aero(n).CLstd = std(CL);
        aero(n).CDstd = std(CD);
        aero(n).CSstd = std(CS);
        aero(n).CRMstd = std(CRM);
        aero(n).CPMstd = std(CPM);
        aero(n).CYMstd = std(CYM);
    end
    %sort by increasing angle of attack
    [~,indx] = sort([aero.AoA]);
    aero = aero(indx);

    processing_file = [mfilename '.m']; %#ok<NASGU>
    processing_date = date; %#ok<NASGU>
    save([src filesep 'Results.mat'],'params','src','aero','tare','processing_file','processing_date');
    Save_TXT(src,params,aero);
end

function Save_TXT(src,params,aero)

    fid = fopen([src,filesep,'Results.txt'],'w+');

    fprintf(fid,'Ref_Area[m2] Ref_Chord_Length[m] Ref_Span_b[m] Xc[m] Yc[m] Zc[m] Balance_AoA_Offset[deg] Model_AoA_Offset[deg] Balance_Side_Slip_Offset[deg] Model_Side_Slip_Slip_Offset[deg] Balance_Roll_Offset[deg] Model_Roll_Offset[deg]\n');
    fprintf(fid,'%s\n',num2str([params.Area,params.Chord,params.span,params.Xmrp,params.Ymrp,params.Zmrp,params.AoAB,params.AoAM,params.AoSB,params.AoSM,params.AoRB,params.AoRM]));

    fprintf(fid,'AoA CL CD CS CRM CPM CYM SD(CL) SD(CD) SD(CS) SD(CRM) SD(CPM) SD(CYM)');

    for n = 1:length(aero)
        fprintf(fid,'\n');  
        fprintf(fid,num2str([aero(n).AoA,aero(n).CL,aero(n).CD,aero(n).CS,aero(n).CRM,aero(n).CPM,aero(n).CYM,aero(n).CLstd,aero(n).CDstd,aero(n).CSstd,aero(n).CRMstd,aero(n).CPMstd,aero(n).CYMstd]));

    end

    fclose(fid);

end

function [cal,f_indx,T_indx] = Calibration(loadcell)
    switch lower(loadcell)
        case 'jr3_fz_down'
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

function [flist,src_dir,folders] = getfiles(condition,varargin)
%[flist,src_dir,folders] = getfiles(condition,varargin)
%Function returns list of files and folders (and source directory) based on
%a given condition.
%Inputs:
%           condition:  conditional string on which returned files are
%                       dependent.  Accepts wildcards.  [] returns all
%                       files in directory.
%           options:
%                       path:   specify top-level directory for searching.
%                               If this is not specified, the program
%                               prompts the user.
%                       '-a' or '-all': searches subdirectories in addition
%                                       to top-level.
%                       '-s' or '-sub': allows user to select a subset of
%                                       the found files.

    if ~exist('condition','var'), condition = ''; end
    if isempty(varargin), varargin = {''}; end

    %Determine search directory
    test = cellfun(@(lst) isdir(lst),varargin);
    if any(test)
        src_dir = varargin{test};
    else
        src_dir = uigetdir('Specify Directory');
    end

    %find files within directory that match given condition
    flist = findfiles(src_dir,condition);

    %find folders and contents within
    cutoff = length(src_dir)+2; %cutoff for folder name
    if any(strcmpi(varargin,'-all')) || any(strcmpi(varargin,'-a'))
        tmp = genpath(src_dir);
        folders = regexp(tmp,pathsep,'split'); %find subdirs
        folders = folders(2:end); %remove top level directory from subdir list
        folder_keep = true(1,length(folders)); %logical check for folder contents
        for n = 1:length(folders)
           tflist = findfiles(folders{n},condition); %find files in given subdir
           if ~isempty(tflist)
               tflist = cellfun(@(file) fullfile(folders{n}(cutoff:end),file),tflist,'UniformOutput',false); %append subdir name to file
               flist = [flist, tflist]; %append new subdir filelist to dir filelist
           else
               folder_keep(n) = false;
           end
        end
        folders = folders(folder_keep);
    else
        folders = cell(0);
    end

    %Lets user choose the subset of files in that directory for gathering
    if any(strcmpi(varargin,'-sub')) || any(strcmpi(varargin,'-s'))
        l = max(cellfun(@length,flist));
        [keep,ok] = listdlg('PromptString','Select Desired Files','SelectionMode','multiple','ListSize',[8*l 300],'ListString',flist);
        if ok==0
            error('Function Terminated Due to User Selection of Cancel')
        end
        flist = flist(keep);
    end
    flist = flist';

    %If user is requesting folders, output file list should be changed
    if nargout > 2
        for n = 1:length(folders)
            group = cellfun(@(x) strcmp(folders{n}(cutoff:end),fileparts(x)),flist);
            swap{n} = flist(group);
        end
        flist = swap;
    end
end

function [tflist] = findfiles(src_dir,condition)
    list = dir([src_dir,filesep,condition]);
    chk = [list.isdir];
    tflist = {list(~chk).name};
end
