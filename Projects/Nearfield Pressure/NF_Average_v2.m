function NF_Average_v2(src_dir,flist,out_dir,savename,cal_dir,params)
%Function acts as an initial processing routine for nearfield/farfield
%acoustic data with trigger signal. Function computes the averaged waveform
%for each channel and then filters this averaged waveform to remove the
%actuator self noise. Outputs are saved in struct labeled 'pf'.
%Code version: 2.0
    
    %Constants
    pf = load(params);
    pf.pp.NCh = length([pf.pp.NFCh pf.pp.FFCh])+1; %total number of channels recorded
    
    %grab file info
    pf.pp.src_dir = src_dir;
    pf.pp.flist = flist;
    pf.pp.cal_dir = cal_dir;
    
    %initialize arrays
    N = length(flist);
    pf.avg_wvfm = cell(1,N);
    pf.sm_wvfm = cell(1,N);
    pf.phys = cell(1,N);    
    
    for n = 1:N  
        disp(['Processing File: ',flist{n}]);
        %create variable name
        [~,fname] = fileparts(flist{n});
        
        %Read filename, get physical parameter data
        pf.phys{n} = NF_ReadFileName_v2(fname,pf.pp);
        
        %Read in data 
        [nfp,ffp,trigger] = ReadData(src_dir,flist{n},pf.pp,pf.phys{n},cal_dir); 
        
        %Average Nearfield waveform
        pf.avg_wvfm{n} = AvgWave([ffp nfp],trigger,pf.pp); %need to fix this line so that the pressure channels do not get reordered if nfp and ffp get reordered
        
        %Filter and smooth averaged nearfield waveforms
        pf.sm_wvfm{n} = NF_wavefilter(pf.avg_wvfm{n},pf.phys{n},pf.pp);
        
    end  
    save([out_dir filesep savename '.mat'],'pf');
end

function [nfp ffp trigger] = ReadData(src_dir,fname,pp,phys,cal_dir)
        %Read file
        fid = fopen([src_dir filesep fname],'r');
        raw = fread(fid,'float32');fclose(fid);
        tmp = reshape(raw,pp.BS,pp.NCh,[]); %reshape data into Points x Channels x Blocks
        
        %Calibrate signals
        tmp = Calibrate(tmp,pp,phys,cal_dir);
        
        %Extract components
        trigger = squeeze(tmp(:,pp.tCh,:)); %extract trigger from rest of voltage traces
        nfp = -tmp(:,pp.NFCh,:); %extract nearfield pressure signal, invert to acount for 4939/2690 mic and conditioner setup
        ffp = -tmp(:,pp.FFCh,:); %extract farfield pressure signal, invert        
end

function [c] = Calibrate(s,pp,phys,cal_dir)
        %Reshape signal for calibration
        %Assumes the signal is in Points x Channels x Blocks
        [N M O] = size(s);
        s = reshape(permute(s,[1 3 2]),N*O,M);

        %Get Microphone Calibrations
        cal = ones(M,1);
        for n = [pp.FFCh pp.NFCh]
            calfile = getfiles(['CAL*Ch' num2str(n) '*mVPa' num2str(cell2mat(phys.gain{n})) '_T*.mat'],cal_dir);
            if length(calfile) > 1, error('Indeterminate Calibration File'); end
            tmp = load([cal_dir filesep calfile{1}]);
            cal(n) = tmp.PaV;
        end
        CAL = spdiags(cal,0,M,M);
        
        %Calibrate signals
        c = s*CAL; %convert pressure from volts to Pa 
        
        %Reshape back into Points x Channels x Blocks
        c = permute(reshape(c,N,O,M),[1 3 2]);
end

function [avg_waveform] = AvgWave(pressure,trigger,pp)
        %Identify actuation indices for each block
        a_indices = cell(1,pp.NB);
        t = zeros(1,pp.NB);
        for n = 1:pp.NB
            a_indices{n} = round(IdentifyActuation(trigger(:,n))); %actuation indices
            t(n) = mean(diff(a_indices{n}));
        end
        T = round(mean(t)); %period of actuation, in indices
        wndo = T; %length of averaging window, in indices

        %Identify actuation windows, sum waveform
        waveform = zeros(wndo,pp.NCh-1,pp.NB);
        for n = 1:pp.NB
            N = sum(a_indices{n}(:)+T < pp.BS); %number of complete waveforms for averaging
            for nn = 1:N
                waveform(:,:,n) = waveform(:,:,n) + pressure(a_indices{n}(nn):a_indices{n}(nn)+wndo-1,:,n)/N;
            end
        end
        avg_waveform = mean(waveform,3);
end

function [phys] = NF_ReadFileName_v2(filename,pp)
%Parses specified filename to determine jet operating conditions (TTR, exit
%velocity, acoustic Mach number, forcing frequency, etc).  Fields are added
%to the 'phys' structure.
    
    md = getMetaFromStr(filename,'NearFieldParams_v2'); %pull metadata from filename
    
    %pull necessary data
    phys.Stdf = md.S.value;
    phys.M = md.M.value;
    phys.T = md.T.value;
    phys.x = md.x.value + (pp.sp-md.x.value)*cosd(md.a.value);%axial distances of microphones 
    phys.y = md.r.value+(phys.x-phys.x(1))*tand(md.a.value); % y/D
    phys.r = zeros(1,length([pp.NFCh pp.FFCh]));
    phys.r(pp.FFCh) = pp.R;
    phys.r(pp.NFCh) = sqrt(phys.x.^2+phys.y.^2)*0.0254; %radial distance of nearfield microphones, converted to meters
    gain_label = [num2str(md.nca.value)';num2str(md.ncb.value)';num2str(md.ncc.value)';num2str(md.ncd.value)';num2str(md.nce.value)'];
    phys.gain = cell(1,pp.NCh-1);
    for n = 1:pp.NCh-1
        phys.gain{n} = GetGain(gain_label(n));
    end
    
    %Calc jet/ambient properties
    phys.To = phys.T + 273.15; %convert to kelvin
    Te = phys.To/(1+1/5*phys.M^2); %jet exit temperature, in kelvin
    phys.TTR = phys.To/(pp.AT+273.15); %jet temperature ratio
    phys.Ue = phys.M*sqrt(1.4*287.05*Te); %nozzle exit velocity (m/s)
    phys.Ma = phys.M*sqrt(Te/(pp.AT+273.15)); %acoustic jet Mach number
    phys.a = sqrt(1.4*287.05*(pp.AT+273.15)); %ambient speed of sound (m/s)
    phys.c = phys.Ue/phys.M; %speed of sound in jet (m/s)
    phys.Uc = phys.Ue*phys.a/(phys.a+phys.c); %estimated convective velocity
    phys.t_c = phys.r/phys.Uc; %convective time for structures
    phys.i_c = round(phys.t_c*pp.FS); %convective index for structures
    phys.t_a = phys.r/phys.a; %acoustic time for actuator pulse
    phys.i_a = round(phys.t_a*pp.FS); %acoustic index for actuator pulse
    end

function [gain] = GetGain(label)
    %gain labels and voltage conversion
    glabel = {'0',  '1',    '2',    '3',    '4',    '5',    '6',    '7',    '8',        '9',        'X'};
    gvalue = {'1',  '3.16', '10',   '31.6', '100',  '316',  '1000', '3160', '10000',    '31600',    '100000'};
    n = strmatch(label,glabel);
    gain = gvalue(n);
end

function [a_i] = IdentifyActuation(trigger)
    %find minimum recorded voltage. This will be assumed to be the drop
    %voltage, as the true voltage is not quite what the user specifies in
    %the waveform generator.
    
    %subtract out mean of trigger to ensure there are zero crossings
    trigger = trigger - mean(trigger);
    vmin = min(trigger);
    
    %calculate first derivative of trigger. This will be used for
    %sub-sample interpolation
    dtrigger = diff(trigger); %first derivative of trigger
    alpha = 0.05*max(dtrigger); %cutoff value for slope determination
    dtrigger(dtrigger < alpha) = []; %get rid of all values not occurring on the rising slope
    slope = mean(dtrigger); %average rise time (per index) of voltage rise
    
    %finds negative peaks in trigger.
    [~,actuation] = findpeaks(-trigger,'minpeakheight',-0.5*vmin,'minpeakdistance',2); %unrefined actuation indices
    actuation = actuation'+double(trigger(actuation-1) > 0); %shifts indices by one if the found peak is on the voltage drop, rather than the voltage rise
    a_i = actuation + (vmin-trigger(actuation))/slope; %interpolate based on average slope of voltage rise
    
    %add fix for incorrect waveform (rather than drop, the trigger rose,
    %then dropped, then returned to zero)
%     vmax = max(trigger);
%     time = round((vmax-vmin)/slope/2); %half width of pulse
%     a_i = a_i-time;
    
    %get rid of any negative indices that may occur due to the
    %interpolation.
    a_i = a_i(a_i > 0.5);
end

function [sm_wave] = NF_wavefilter(avg_wvfm,phys,pp)
    %Filters and smooths the averaged nearfield waveforms in order to
    %remove the actuator self-noise

    [N,M] = size(avg_wvfm);

    %wavelet params
    mother = 'paul'; %mother wavelet
    param = 4; %wavelet parameter, mother-specific
    dt = 1/pp.FS; %sampling time
    s0 = 1/pp.FS; %smallest scale
    J1 = 4e2; %number of scales (minus 1)
    DJ = log2(3*N/2)/J1; %scale spacing, J1 = (LOG2(N DT/S0))/DJ.
    
    %smoothing params
    wt = 20; %search window temporal width 
    hww = 15; %hard wavelet smoothing window
    sww = 7; %soft wavelet smoothing window
    htw = 15; %hard temporal smoothing window
    stw = 5; %soft temporal smoothing window
    hwf = @hamming;
    swf = @hann;
    htf = @hamming;
    stf = @hann;
    
    %shift signal so acoustic time of arrival is centered in window
    center = ceil(N/2);
    i_a = rem(phys.i_a,N); %remove periods when acoustic indice > window length 
    shift = center - i_a;
    
    %initialize output matrix
    sm_wave = zeros(size(avg_wvfm));
    
    for n = 1:M
        %rotate signal so expected time of arrive is in center of window
        s = avg_wvfm(:,n);
        s = circshift(s,shift(n));
        
        %normalize signal
        sm = mean(s); %mean for channel
        sn = s-sm; %normalized signal
           
        %wavelet transform
        [wave,~,scale,~,~, ~, k] = contwt(sn,dt,0,DJ,s0,J1,mother,param);
        dims = size(wave);
        
        %find actuator self-noise window
        tw = [center-wt center+wt]; %temporal index of window 
        if tw(1) < 1, tw(1) = 2; end
        if tw(2) > dims(2), tw(2) = dims(2)-1; end
        pks = imregionalmax(abs(wave(:,tw(1)-1:tw(2)+1))); %find local maxima
        pks = pks(:,2:end-2); %get rid of values on edge of window (they are not true peaks)
        [ps pt] = find(pks); %get temporal and scale indices for each of the peaks
        pt = pt + tw(1)-1; %shift window so it matches original location
        pvals = diag(abs(wave(ps,pt))); %peak values
        pmean = mean(abs(wave(ps,:)),2); %peak mean
        pstd = std(abs(wave(ps,:)),1,2); %peak std
        if ~isempty(pvals)
            chk = pvals > pmean + pstd; %check to make sure peaks are statistically important
            pt = pt(chk); %throw out unimportant peaks
            ps = ps(chk);
            ias = [1 max(ps)+floor(wt/2)]; %scale window size
            iat = [min(pt)-floor(wt/2) max(pt)+floor(wt/2)]; %temporal window size  
        end

        %filter out actuator self noise in wavelet domain
        fwave = wave; %filtering wavelet coefs
        if ~isempty(ps)
            if iat(1) <= hww, iat(1) = hww+1; end
            if dims(2)-iat(2) <= hww, iat(2) = dims(2)-hww; end 
            hws = ndnanfilter(wave(ias(1):ias(2)+hww,iat(1)-hww:iat(2)+hww),hwf,[hww hww]); %hard smoothing in wavelet domain
            fwave(ias(1):ias(2),iat(1):iat(2)) = hws(1:ias(2),hww+1:end-hww); %replace window
            fwave = ndnanfilter(fwave,swf,[sww sww],[],{},{'circular'}); %soft smoothing in wavelet domain
        end
                
        %reconstruct waveform
        filtered = s;
        if ~isempty(ps)            
            rec = invcwt(fwave,mother, scale,param,k)+sm; %reconstruct waveform from filtered wavelet coefs, add in signal mean
            filtered(iat(1):iat(2)) = rec(iat(1):iat(2)); %replace window
        end
        
        %smooth final waveform in temporal domain
        if ~isempty(ps)
            if iat(1) <= htw, iat(1) = htw+1; end
            hts = wmavg(filtered(iat(1)-htw:iat(2)+htw),htw,htf,1); %hard smoothing in temporal domain
            filtered(iat(1):iat(2)) = hts(htw+1:end-htw); %replace window
        end
        sm = wmavg(filtered,stw,stf,1);
        
        %un-shift window
        sm_wave(:,n) = circshift(sm,-shift(n));
    end  
end