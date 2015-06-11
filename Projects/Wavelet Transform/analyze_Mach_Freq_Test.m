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

%   Specify the base directory to be searched for data files
base_dir = 'G:\AniLabData\Expt_Pr\20110824-31 - MachFreqSweep\20110830_MachFreqSweep\1inch_repeat';

%   Jet measurement conditions - specified
pressure_loc_x = (2:0.5:4)*0.0254;   %Axial locations of pressure transducers
t_sample = 5e-6;            %Sampling period

%   Jet measurement conditions - derived
n_pressure = length(pressure_loc_x);%# pressure signals recorded

%   The data for different Mach numbers are stored in folders with names of
%   the form 'M0.9' say, for Mach 0.9.
temp_files = dir(fullfile(base_dir,'M*'));  %Read all files/folders starting with 'M'
temp_files_n = length(temp_files);
Machs = zeros(temp_files_n,1);              %Pre-allocate
Mach_strs = cell(temp_files_n,1);           %Pre-allocate
n_Machs = 0;                                %# of valid folders
for mdx = 1:temp_files_n
    fn = temp_files(mdx).name;
    if isdir(fullfile(base_dir,fn)) && ~isnan(str2double(fn(2:end)))
        %It's a folder with expected naming
        n_Machs = n_Machs + 1;
        Machs(n_Machs) = str2double(fn(2:end));
        Mach_strs{n_Machs} = fn(2:end);
    end
end
Machs = Machs(1:n_Machs);
[Machs,Mach_srt_idx] = sort(Machs);
Mach_strs = Mach_strs(Mach_srt_idx);

%   Obtain # freqs tested from processed data for the first Mach No.
Mach_dir1 = fullfile(base_dir,['M',Mach_strs{1}]);
temp = load(fullfile(Mach_dir1,['M',Mach_strs{1},'_avg_waveform.mat']));
n_freqs = length(temp.F_f);

wvfrm_params.F_f = temp.F_f;
wvfrm_params.F_f_str = temp.F_f_str;
wvfrm_params.Machs = Machs;
wvfrm_params.Mach_strs = Mach_strs;

clear temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Determine the parameters of the average pressure waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wvfrm_params.time_2lo = zeros(n_Machs,n_pressure,n_freqs);
wvfrm_params.time_lo2hi = zeros(size(wvfrm_params.time_2lo));
wvfrm_params.amp_hi = zeros(size(wvfrm_params.time_2lo));
wvfrm_params.amp_lo = zeros(size(wvfrm_params.time_2lo));
wvfrm_params.U_j = zeros(n_Machs,1);

mother = 'PAUL';
order = 4;
FS = 1/t_sample;
NS = 400; %number of scales
NS = NS-1;
ss = 1/FS;	%smallest scale - sec
LS = 1e-2;	%Largest scale - sec
SS = log2(LS/ss)/NS;	%Scale spacing - sec

for mdx = 1:n_Machs
    Mach_dir = fullfile(base_dir,['M',Mach_strs{mdx}]);
    temp = load(fullfile(Mach_dir,['M',Mach_strs{mdx},'_avg_waveform.mat']));
    avg_waveform = temp.avg_waveform;
    t_waveform = temp.t_waveform;
    F_f = temp.F_f;
    wvfrm_params.U_j(mdx) = mean(temp.U_j);
    clear temp
    if length(F_f) ~= n_freqs
        disp(['# F_f''s doesn''t match for M = ',Mach_strs{mdx}]);
        continue;
    end
    for fdx = 1:n_freqs
        T_f = 1/F_f(fdx);
        if T_f < 1.5*diff(t_waveform([1,end],1)) && fdx > 1
            n_sample_T_f = ceil(T_f/t_sample);
            t_wvfrm_curr = zeros(n_sample_T_f,n_pressure);
            for pdx = 1:n_pressure
                mean_time_2lo = mean(wvfrm_params.time_2lo(mdx,pdx,1:fdx-1));
                [~,index_time_2lo] = min(abs(t_waveform(:,pdx) - mean_time_2lo));
                index_first = index_time_2lo - round(n_sample_T_f/2);
                index_last = index_first + 2*n_sample_T_f - 1;
                avg_waveform(1:index_first-1,pdx,fdx) = 0;
                avg_waveform(index_last+1:end,pdx,fdx) = 0;
            end
        end
        for pdx = 1:n_pressure
            [wvfrm_params.amp_lo(mdx,pdx,fdx),idx_lo] = min(avg_waveform(:,pdx,fdx));
            [wvfrm_params.amp_hi(mdx,pdx,fdx),idx_hi] = max(avg_waveform(:,pdx,fdx));
            wvfrm_params.time_2lo(mdx,pdx,fdx) = t_waveform(idx_lo,pdx);
%             wvfrm_params.time_lo2hi(mdx,pdx,fdx) = diff(t_waveform([idx_lo,idx_hi],pdx));
%             if wvfrm_params.time_lo2hi(mdx,pdx,fdx) < 0
%                 wvfrm_params.time_lo2hi(mdx,pdx,fdx) = T_f + wvfrm_params.time_lo2hi(mdx,pdx,fdx);
%             end
            [WAVE,PERIOD,SCALE,COI,tmp,tmp2,k] = contwt(avg_waveform(:,pdx,fdx),1/FS,0,SS,ss,NS,mother,order);
            [~,scale_idx] = max(max(imag(WAVE).^2,[],2),[],1);
            wvfrm_params.time_lo2hi(mdx,pdx,fdx) = PERIOD(scale_idx);
        end
    end
end

save(fullfile(base_dir,'avg_waveform_parameters1.mat'),'wvfrm_params');
