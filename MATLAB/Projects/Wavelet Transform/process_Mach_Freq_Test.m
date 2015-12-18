%   Analyzes the results from INPUT/OUTPUT tests with varying frequencies
%   and Mach numbers.

%   These tests are conducted using LabVIEW running on 2 computers.
%   The near-field pressure measurements were acquired using one LabVIEW
%   program running on a computer. It also acquired a ramp signal that is
%   synchronized with the actuator control pulse train. The latter is
%   generated in a 2nd computer.

%   The near-field measurements consist of 5 mics on a uniform linear array
%   with 0.5 inches between them.

%   All the 8 actuators are fired simultaneously.

%   A separate program is called to interpret the ramp signal to determine
%   the time instants in the data record when the LAFPA was fired. In
%   general, the firing instant would not coincide with a sampling instant.

%   The goal here is to find a conditional-averaged waveform. The condition
%   is the beginning of actuator firing. Relative to this t=0 instant, the
%   program determines the average pressure signal at the various
%   microphones. 

clear all

%   Load microphone calibration file
Calib_dir = 'G:\AniLabData\Expt_Pr\20110824-31 - MachFreqSweep\20110913_Calibration';
Calib_fn = 'Calib_3mVPa';

%   Specify the base directory to be searched for data files
base_dir = 'G:\AniLabData\Expt_Pr\20110824-31 - MachFreqSweep\20110830_MachFreqSweep\1inch_repeat\M0.9';

temp = load(fullfile(Calib_dir,Calib_fn));
Calib = temp.Calib; %Pa/V per channel

%   Specify Mach No. (this is used in selecting appropriate data files)
Mach_e_str = '0.9'; %String (for data filename matching during retrieval)
Mach_e = str2double(Mach_e_str);   %Jet exit Mach number

%   Jet operating conditions - specified
D_e = 0.0254;   %Jet exit diameter
T_amb = 295;    %Ambient temperature (K)
gas_const_air = 287; %Gas constant for air (J/kg K)
gamma = 1.4;    %Ratio of specific heats for air
TTR = 1;        %Total temperature ratio

%   Jet measurement conditions - specified
pressure_loc_x = (2:0.5:4)*0.0254;   %Axial locations of pressure transducers
t_sample = 5e-6;            %Sampling period
blk_sz = 81920;             %Size of block of contiguous samples
n_blks = 10;                %# contiguous blocks
n_far_field = 7;            %# far-field signals preceding near-field signals

%   Actuation signal to ramp conversion paramters - specified
act_sig_param.sampling_rate = 1/t_sample;
act_sig_param.actuation_ramp.Vpp = 5;
act_sig_param.actuation_ramp.clk_rate = 4e6;
act_sig_param.actuation_ramp.period = 5.0000e-005;

%   Jet operating conditions - derived
a_amb = sqrt(gamma*gas_const_air*T_amb);    %Ambient speed of sound
T_0 = TTR*T_amb;                            %Stagnation temperature
T_e = T_0/(1+(gamma-1)/2*Mach_e^2);         %Exit temperature
a_e = sqrt(gamma*gas_const_air*T_e);        %Exit speed of sound
U_e = a_e*Mach_e;                           %Exit velocity
U_c = U_e*(a_amb/(a_amb+a_e));              %Theoretical convective velocity

%   Jet measurement conditions - derived
n_pressure = length(pressure_loc_x);%# pressure signals recorded
n_ch_save = n_far_field+n_pressure+1;	%# channels saved (actuation signal also)
%Theoretical convective times to different pressure transducers. This
%determines center time for correlation.
t_conv = pressure_loc_x/U_c;
%Expected duration of pressure pulse (be very liberal). This determines
%time window for correlation.
t_event = 15*(D_e/U_c);
n_event = floor(t_event/t_sample);  %# samples in correlation window
%   For ease of later handling, covert n_event to the nearest odd number.
if n_event == 2*round(n_event/2)
    n_event = n_event + 1;
end
max_blks_per_read = 1; %Too many blocks can't be read at once (low memory)
n_reads = ceil(n_blks/max_blks_per_read); %# reads to get all blocks

%   The microphones are located at different axial locations, so that the
%   convective delay from the nozzle exit to the microphone tip is
%   different. 
%   'window_indices' is a 2-D array of indices numbered from the instant of
%   any pulse firing. Its columns are for the different microphones. Its
%   rows have 2 entries: the start and end indices from the beginning to
%   the end of the correlation window for the particular microphone.
%   For each microphone, the window is centered on the index corresponding
%   to the expected convective delay from the nozzle exit to the particular
%   mic.
window_indices = zeros(2,n_pressure);
t_waveform = zeros(n_event,n_pressure);
for pdx = 1:n_pressure
    N_conv_pdx = round(t_conv(pdx)/t_sample);
    window_indices(:,pdx) = [-(n_event-1)/2,(n_event-1)/2] + N_conv_pdx;
    t_waveform(:,pdx) = (window_indices(1,pdx):window_indices(end,pdx))*t_sample;
end

%   Find all files for forced data
all_files = dir(fullfile(base_dir,['M',Mach_e_str,'*.bin']));
n_files_prelim = length(all_files);
index_set = zeros(n_files_prelim,1);
n_forced_files = 0;
for idx = 1:n_files_prelim
    info = parse_pressure_fname(all_files(idx).name);
    if ~info.baseline
        n_forced_files = n_forced_files + 1;
        index_set(n_forced_files) = idx;
    end
end
forced_files = all_files(index_set(1:n_forced_files));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Determine the average pressure waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg_waveform = zeros(n_event,n_pressure,n_forced_files); %Pre-allocate
F_f = zeros(n_forced_files,1);
F_f_str = cell(n_forced_files,1);
U_j = zeros(n_forced_files,1);
out_fn = fullfile(base_dir,['M',Mach_e_str,'_avg_waveform.mat']);
for fdx = 1:n_forced_files %For each data file stored
%     disp(['Processing file ',num2str(fdx),' of ',num2str(n_forced_files)]);
    n_events_pdx = zeros(n_pressure,1); %# events recorded by individual mics
    info = parse_pressure_fname(forced_files(fdx).name);
    U_j(fdx) = info.thermo.U_j;
    fid = fopen(fullfile(base_dir,forced_files(fdx).name),'r','ieee-le');
    for rdx = 1:n_reads
        n_blks2read_curr = min(n_blks-(rdx-1)*max_blks_per_read,max_blks_per_read);
        dat = fread(fid,blk_sz*n_ch_save*n_blks2read_curr,'float32');
        dat = reshape(dat,blk_sz,n_ch_save,[]);
        dat = permute(dat,[1,3,2]);
        sig_act = dat(:,:,end);
        dat = dat(:,:,n_far_field+(1:n_pressure));
        [~,st_ids,F_f_curr] = Continuous_actuation_signal_decode(sig_act,act_sig_param);
        for bdx = 1:n_blks2read_curr
            st_id_curr_blk = st_ids{bdx};
            for pdx = 1:n_pressure  %For each mic
                %   Only use those events which afford a full window of data
                %   This involves ignoring those events whose corresponding windows
                %   would have the beginning or ending truncated.
                st_id_curr_blk_pdx = st_id_curr_blk(floor(st_id_curr_blk) > - window_indices(1,pdx) ...
                    & ceil(st_id_curr_blk) <= blk_sz - window_indices(end,pdx));
                n_events_curr = length(st_id_curr_blk_pdx);   %# events considered
                n_events_pdx(pdx) = n_events_pdx(pdx) + n_events_curr; %Running count
                for edx = 1:n_events_curr
                    event_beg = st_id_curr_blk_pdx(edx) + window_indices(1,pdx);
                    event_end = st_id_curr_blk_pdx(edx) + window_indices(end,pdx);
                    orig_idx = floor(event_beg):ceil(event_end);
                    avg_waveform(:,pdx,fdx) = avg_waveform(:,pdx,fdx) ...
                        + interp1(orig_idx,dat(orig_idx,bdx,pdx),event_beg:event_end)';
                end
            end
        end
        clear sig_act dat
    end
    fclose(fid);
    F_f(fdx) = F_f_curr;
    F_f_str{fdx} = num2str(F_f_curr);
    for pdx = 1:n_pressure
        avg_waveform(:,pdx,fdx) = avg_waveform(:,pdx,fdx)*(Calib(pdx)/n_events_pdx(pdx));
    end
    save(out_fn,'avg_waveform','t_waveform','F_f','F_f_str','U_j','fdx');
end

%   Sort by forcing frequency
[F_f,F_f_srt_idx] = sort(F_f);
F_f_str = F_f_str(F_f_srt_idx);
U_j = U_j(F_f_srt_idx);
avg_waveform = avg_waveform(:,:,F_f_srt_idx);

%   Determine the legend labels
legend_labels = cell(n_pressure,1);
for pdx = 1:n_pressure
    legend_labels{pdx} = ['x/D = ',num2str(pressure_loc_x(pdx)/D_e,'%.1f')];
end

save(out_fn,'avg_waveform','t_waveform','legend_labels','F_f','F_f_str','U_j');
