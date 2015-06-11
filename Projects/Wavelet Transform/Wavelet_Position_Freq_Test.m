%   Analyzes the results from INPUT/OUTPUT tests with varying frequencies
%   and Mach numbers

%   The goal is to determine various parameters of the phase-averaged
%   impulse response

clear all

%   Specify the base directory to be searched for data files
base_dir = 'C:\Users\Hind Alkandry\Documents\Data\Near Field\Linear traverse (8 Mics)\20111030_M0.9\Raw Data\M0.9\Forced Files\Amp_3';

%   Specify threshold factor of wavelet energy to reject before
%   reconstruction. Setting to 0 means no wavelet filtering.
wvlt_rjct_fact = 0.25;

%   The data for different Mach numbers are stored in folders with names of
%   the form 'M0.9' say, for Mach 0.9.
temp_files = dir(fullfile(base_dir,'P*'));  %Read all files/folders starting with 'M'
temp_files_n = length(temp_files);
Posi = zeros(temp_files_n,1);              %Pre-allocate
Posi_strs = cell(temp_files_n,1);          %Pre-allocate
n_Posi = 0;                                %# of valid folders
for mdx = 1:temp_files_n
    fn = temp_files(mdx).name;
    if isdir(fullfile(base_dir,fn)) && ~isnan(str2double(fn(2:end)))
        %It's a folder with expected naming
        n_Posi = n_Posi + 1;
        Posi(n_Posi) = str2double(fn(2:end));
        Posi_strs{n_Posi} = fn(2:end);
    end
end
Posi = Posi(1:n_Posi);
[Posi,Posi_srt_idx] = sort(Posi);
Posi_strs = Posi_strs(Posi_srt_idx);

wvfrm_params.Posi = Posi;
wvfrm_params.Posi_strs = Posi_strs;

%   Obtain # freqs tested from processed data for the first Mach No.
Posi_dir1 = fullfile(base_dir,['P',Posi_strs{1}]);
temp = load(fullfile(Posi_dir1,['P',Posi_strs{1},'_avg_waveform.mat']));
avg_wvfrm = temp.avg_wvfrm;
clear temp
n_freqs = length(avg_wvfrm.f_F);
n_pressure = length(avg_wvfrm.pressure_loc_x);
t_sample = diff(avg_wvfrm.t_waveform(1:2,1));

wvfrm_params.f_F = avg_wvfrm.f_F;
wvfrm_params.f_F_str = avg_wvfrm.f_F_str;

clear avg_wvfrm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reconstruct the average pressure waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wvfrm_params.U_j = zeros(n_Posi,1);

mother = 'PAUL';
order = 4;
FS = 1/t_sample;
NS = 400; %number of scales
NS = NS-1;
ss = 1/FS;	%smallest scale - sec
LS = 1e-2;	%Largest scale - sec
SS = log2(LS/ss)/NS;	%Scale spacing - sec

for mdx = 1:n_Posi
    Posi_dir = fullfile(base_dir,['P',Posi_strs{mdx}]);
    temp = load(fullfile(Posi_dir,['P',Posi_strs{mdx},'_avg_waveform.mat']));
    avg_wvfrm = temp.avg_wvfrm;
    clear temp
    
    avg_waveform = avg_wvfrm.avg_waveform;
    t_waveform = avg_wvfrm.t_waveform;
    f_F = avg_wvfrm.f_F;
    wvfrm_params.U_j(mdx) = mean(avg_wvfrm.U_j);
    clear avg_wvfrm
    
    n_samples = size(avg_waveform,1);
    if length(f_F) ~= n_freqs
        disp(['# f_F''s doesn''t match for P= ',Posi_strs{mdx}]);
        continue;
    end
    for fdx = 1:n_freqs
        T_f = 1/f_F(fdx);
        
        %wvfrm_recon = zeros(size(avg_waveform));
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
                wvfrm_params.wvfrm_recon(:,pdx) = invcwt(WAVE,mother,SCALE,order,k);
            else
                wvfrm_params.wvfrm_recon(:,pdx) = wvfrm_curr;
            end
                wvfrm_params.avg_waveform= avg_waveform;
                wvfrm_params.t_waveform= t_waveform;
        end
    end


%   Specify the base directory for the data to be saved to
save_dir = 'C:\Users\Hind Alkandry\Documents\Data\Near Field\Linear traverse (8 Mics)\20111030_M0.9\Processed Data\20120308 (Ani Code)';
out_fn = fullfile(save_dir,['P',Posi_strs{mdx},'avg_wvfrm_params_',num2str(wvlt_rjct_fact),'.mat']);
save(out_fn,'wvfrm_params');
end