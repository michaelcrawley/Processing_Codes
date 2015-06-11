%   Analyzes the results from INPUT/OUTPUT tests with varying frequencies
%   and Mach numbers

%   The goal is to determine various parameters of the phase-averaged
%   impulse response

clear all

%   Specify the base directory to be searched for data files
base_dir = 'C:\Users\Aniruddha\Documents\GDTL_Docs\AniLabWork\20110830_MachFreqSweep\1.5inch';

%   Specify threshold factor of wavelet energy to reject before
%   reconstruction. Setting to 0 means no wavelet filtering.
wvlt_rjct_fact = 0.25;

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

wvfrm_params.Machs = Machs;
wvfrm_params.Mach_strs = Mach_strs;

%   Obtain # freqs tested from processed data for the first Mach No.
Mach_dir1 = fullfile(base_dir,['M',Mach_strs{1}]);
temp = load(fullfile(Mach_dir1,['M',Mach_strs{1},'_avg_waveform.mat']));
avg_wvfrm = temp.avg_wvfrm;
clear temp
n_freqs = length(avg_wvfrm.f_F);
n_pressure = length(avg_wvfrm.pressure_loc_x);
t_sample = diff(avg_wvfrm.t_waveform(1:2,1));

wvfrm_params.f_F = avg_wvfrm.f_F;
wvfrm_params.f_F_str = avg_wvfrm.f_F_str;

clear avg_wvfrm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Determine the parameters of the average pressure waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wvfrm_params.time_2lo = zeros(n_Machs,n_pressure,n_freqs);
wvfrm_params.time_lo2hi = zeros(size(wvfrm_params.time_2lo));
wvfrm_params.time_lo20 = zeros(n_Machs,n_pressure,n_freqs);
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
    avg_wvfrm = temp.avg_wvfrm;
    clear temp
    
    avg_waveform = avg_wvfrm.avg_waveform;
    t_waveform = avg_wvfrm.t_waveform;
    f_F = avg_wvfrm.f_F;
    wvfrm_params.U_j(mdx) = mean(avg_wvfrm.U_j);
    clear avg_wvfrm
    
    n_samples = size(avg_waveform,1);
    if length(f_F) ~= n_freqs
        disp(['# f_F''s doesn''t match for M = ',Mach_strs{mdx}]);
        continue;
    end
    for fdx = 1:n_freqs
        T_f = 1/f_F(fdx);

        for pdx = 1:n_pressure
            wvfrm_curr = avg_waveform(:,pdx,fdx);

            %   Perform wavelet filtering if required
            if wvlt_rjct_fact > eps
                %   Transform to wavelet domain
                [WAVE,PERIOD,SCALE,COI,tmp,tmp2,k] = contwt(wvfrm_curr,1/FS,0,SS,ss,NS,mother,order);

                %   Filter in the wavelet domain
                ener = abs(WAVE).^2;
                max_ener = max(max(ener));
                WAVE(ener < max_ener*wvlt_rjct_fact) = 0;

                %   Reconstruct the response from filtered wavelet domain
                wvfrm_recon = invcwt(WAVE,mother,SCALE,order,k);
            else
                wvfrm_recon = wvfrm_curr;
            end

            %   If more than one forcing period is captured in the
            %   averaging window (which cannot be allowed for the first
            %   (lowest) forcing frequency case), then make the
            %   phase-averaged waveform go to 0 outside one forcing period.
            %   The temporal location of the phase-averaged waveform to
            %   retain is determined such that the average location of the
            %   negative peaks from all previous forcing frequency cases
            %   falls at the quarter mark of the retained window. This is
            %   done from observation of the waveforms: they line up at
            %   about the negative peak for most forcing frequencies, and
            %   for continuously oscillating responses (i.e. no quiet
            %   period) the negative peak falls at about the quarter mark.
            if diff(t_waveform([1,end],1)) > T_f && fdx > 1
                %# samples in one forcing period
                n_sample_T_f = ceil(T_f/t_sample);
                %Average time-of-arrival of -ve peak from all lower forcing
                %frequencies
                mean_time_2lo = mean(wvfrm_params.time_2lo(mdx,pdx,1:fdx-1));
                %Index corresponding to this time
                [~,index_time_2lo] = min(abs(t_waveform(:,pdx) - mean_time_2lo));
                %Index of last instant of response to set to 0 starting
                %from the first
                index_first = max(index_time_2lo - round(n_sample_T_f/4),1);
                %Index of last instant of response to retain as is with all
                %instants set to 0 beyond this
                index_last = min(index_first + round(n_sample_T_f) - 1,n_samples);
                %Do the actual setting to 0
                wvfrm_recon(1:index_first-1) = 0;
                wvfrm_recon(index_last+1:end) = 0;
            end
            
            %   Find index of instant for minimum (should be -ve peak)
            [~,t_idx_lo] = min(wvfrm_recon);
            %   Record corresponding time as time-of-arrival
            wvfrm_params.time_2lo(mdx,pdx,fdx) = t_waveform(t_idx_lo,pdx);
            %   Record corresponding amplitude as minimum
            wvfrm_params.amp_lo(mdx,pdx,fdx) = wvfrm_curr(t_idx_lo);
            
            %   Set response at all instants prior to -ve peak as zero.
            %   This ensures that the +ve peak will be picked from AFTER
            %   the -ve peak
            wvfrm_recon(1:t_idx_lo) = 0;
            
            %   To determine 0-crossing
            t_idx_0 = find(wvfrm_recon > 0,1,'first');
            if ~isempty(t_idx_0)
                %Interpolate between two points on either side of 0 to
                %determine precise 0-crossing time
                wvfrm_params.time_lo20(mdx,pdx,fdx) ...
                    = interp1(wvfrm_recon(t_idx_0+(-1:0)), ...
                        t_waveform(t_idx_0+(-1:0),pdx),0) ...
                        - t_waveform(t_idx_lo,pdx);
            end
            
            %   Retain only +ve amplitudes, and default the time_lo2hi to 0
            %   in case no +ve amplitude remains. This avoids spurious
            %   cases where the -ve peak appears after the +ve peak
            wvfrm_recon(wvfrm_recon < 0) = 0;
            if norm(wvfrm_recon) < eps
                continue;
            end
            
            %   Find index of instant for maximum
            [~,t_idx_hi] = max(wvfrm_recon);
            %   Determine time from -ve peak to +ve peak
            wvfrm_params.time_lo2hi(mdx,pdx,fdx) = diff(t_waveform([t_idx_lo,t_idx_hi],pdx));
            %   Determine amplitude of +ve peak
            wvfrm_params.amp_hi(mdx,pdx,fdx) = wvfrm_curr(t_idx_hi);
        end
    end
end

save(fullfile(base_dir,['avg_waveform_parameters_',num2str(wvlt_rjct_fact),'.mat']),'wvfrm_params');
