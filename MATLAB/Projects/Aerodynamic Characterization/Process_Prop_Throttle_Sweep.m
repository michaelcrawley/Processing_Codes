function loads = Process_Prop_Throttle_Sweep(src,data_files,tare_file,loadcell)

    [cal,f_indx] = LoadCell_Calibration(loadcell);
    
    %Process tare
    voltage = Read_LoadCell([src,filesep,tare_file],9);
    tare_forces = (cal*voltage)';
    tare_mean = mean(tare_forces);
    tare_std = std(tare_forces);
    
    %Process data files
    num_data = length(data_files);
    loads(1:num_data) = struct();
    for n = 1:num_data
        %Grab throttle percentage & tunnel velocity
        conditions = regexp(data_files{n},'_WT(?<vel>[0-9.-]+)+_THR(?<throttle>[0-9.-]+)','names');
        loads(n).velocity = str2double(conditions.vel);
        loads(n).throttle = str2double(conditions.throttle);
        
        %Read raw voltage and convert to loads
        voltage = Read_LoadCell([src,filesep,data_files{n}],9);
        raw = (cal*voltage)';
        
        %Remove tare, calc SD & SE convert to thrust
        corrected = raw - repmat(tare_mean,size(raw,1),1);
        forces = f_indx(corrected);
                
        %%%%%%%%%%%%Hard-coding orientation of loadcell for now
        downsample = 400; % natural frequency of the setup is ~5 Hz, and sample rate was 2000 Hz
        num_independent_samples = length(forces(1:downsample:end,1));
        SD = sqrt(std(forces).^2 + f_indx(tare_std).^2); %standard deviation along load cell axes
        %Calc mean thrust
        forces = mean(forces);
        loads(n).thrust = sqrt(sum(forces(1:2).^2));
        loads(n).axis = atand(forces(2)/forces(1));
        %Calc thrust standard error
        thrust_SD = 0.5*(forces(1)^2 + forces(2)^2)^(-1/2)*sqrt((2*forces(1)*SD(1))^2 + (2*forces(2)*SD(2))^2);
        loads(n).SE = thrust_SD/sqrt(num_independent_samples);        
    end
end