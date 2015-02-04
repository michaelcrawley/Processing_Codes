function NF_Preprocessing(src_dir,flist,cal_dir,out_dir,params,cpus)
    %This is a test function, the purpose of which is to preprocess the 
    %data files. For the forced cases, it will filter out the
    %actuator self noise of the nearfield microphone signals before
    %phase-averaging (as opposed to the previous method where it was done
    %afterwards). 
    %Code version: 1.1
    %Output mat-file version: 1.1
    
    %Constants
    load(params);
%     pp.BS = 81920; %block size of data
%     pp.sBS = 8192; %sub-block size
%     pp.NB = 10;
%     
%     pp.sNBmax = 100;
%     pp.NFCh = 13:20; %nearfield microphone channels
%     pp.lNFCh = 12:19; %nearfield channels after trigger has been removed
%     pp.FFCh = 1:11; %farfield microphone channels    
%     pp.tCh = 12; %channel number of trigger signal
%     
%     pp.FS = 200000; %data acquisition frequency, in Hz
%     pp.D = 0.0254;
%     pp.AT = 26.1; %atmospheric temperature, in C
%     pp.sp = 1;
%     pp.nm = 8;
%     pp.R = [2.45 2.57 2.69 2.90 3.19 3.46 3.68 3.28 3.15 2.63 2.56]; %radial distance of farfield microphones, in meters
    pp.sNBmax = floor(pp.BS/pp.sBS)*pp.NB;
    pp.NCh = length([pp.NFCh pp.FFCh])+1; %total number of channels recorded
    %make output directory
    mkdir(out_dir);
    
    %open workers
    if ~exist('cpus','var'), cpus = 1; end
    if matlabpool('size') > 0
        matlabpool close;
    end
    matlabpool(cpus);
    
    %initialize variable
    tl = cell(1,length(flist));

    parfor n = 1:length(flist)
        disp(['Processing File: ',flist{n}]);
        %create variable name
        filename = flist{n};
        [~,fname] = fileparts(filename);        
        
        %Read filename, get physical parameter data
        tl{n}.phys = NF_ReadFileName_v2(fname,pp);
        
        %grab processing info
        tl{n}.src_dir = src_dir;
        tl{n}.filename = filename;
        
        %Read in data 
        [nfp,ffp,tl{n}.trigger.sig] = ReadData(src_dir,flist{n},pp,tl{n}.phys,cal_dir); 
                        
        %reshape into processing sub-blocks
        tl{n}.nf.pblocks.p = permute(nfp, [1 3 2]);
        tl{n}.ff.pblocks.p = permute(ffp, [1 3 2]);
        tl{n}.nf.pblocks.p = reshape(tl{n}.nf.pblocks.p,pp.sBS,[],length(pp.NFCh));
        tl{n}.ff.pblocks.p = reshape(tl{n}.ff.pblocks.p,pp.sBS,[],length(pp.FFCh));
        tl{n}.nf.pblocks.p = permute(tl{n}.nf.pblocks.p, [ 1 3 2]);        
        tl{n}.ff.pblocks.p = permute(tl{n}.ff.pblocks.p, [ 1 3 2]);
        
        if tl{n}.phys.Stdf ~= 0  %forced cases
            %Reshape data into forcing sub-blocks
            [pressure tl{n}.trigger.a_i]= iSubBlocks([ffp nfp],tl{n}.trigger.sig,pp);           

            %Smooth each individual sub-block
            [~,~,NSB] = size(pressure);
            sm_pressure = zeros(size(pressure));
            for nn = 1:NSB                
                %print processing completeness
                fprintf('\t %d of %d \n',[nn NSB]);
                                
                %smooth block
                sm_pressure(:,:,nn) = wavefilter(pressure(:,:,nn),tl{n}.phys,pp);
            end

            %Reshape raw/filtered signal into forcing sublocks
            tl{n}.nf.ablocks.p = pressure(:,pp.lNFCh,:);
            tl{n}.ff.ablocks.p = pressure(:,pp.FFCh,:);
            tl{n}.nf.ablocks.smp = sm_pressure(:,pp.lNFCh,:);
            tl{n}.ff.ablocks.smp = sm_pressure(:,pp.FFCh,:);
            
            %Test new conversion function
            tl{n}.test = Convert_AtP(tl{n}.trigger.a_i,[ffp nfp],sm_pressure,pp);

            %break into processing sub-blocks           
            slocator = 1;
            locator = 1;
            sm_pressure = permute(sm_pressure,[1 3 2]);
            smp = zeros(pp.sBS,pp.sNBmax,pp.NCh-1);
            for nn = 1:pp.NB
                nsblocks = length(tl{n}.trigger.a_i{nn});
                tmp = sm_pressure(:,slocator:(nsblocks+slocator-1),:);
                [X,Y,Z] = size(tmp);
                gb = rem(X*Y,pp.sBS);
                tmp = reshape(tmp,[X*Y Z]);
                tmp = tmp(1:end-gb,:);
                nblocks = floor(length(tmp)/pp.sBS);
                tmp = reshape(tmp,[pp.sBS nblocks Z]);
                smp(:,locator:(nblocks+locator-1),:)=tmp;
                
                locator = locator+nblocks;
                slocator = slocator+nsblocks;
            end
            smp = permute(smp,[1 3 2]);
            smp = smp(:,:,1:(locator-1)); %remove blank blocks
            tl{n}.nf.pblocks.smp = smp(:,pp.lNFCh,:);
            tl{n}.ff.pblocks.smp = smp(:,pp.FFCh,:);
            sm_pressure = []; smp = [];           

        end 
        SaveData(tl{n},[out_dir filesep fname '_proc.mat']);
        tl{n} = [];
    end
    matlabpool('close');
end

function [output] = Convert_AtP(a_i,p_p,sm_p_a,pp)
    %This function will reorder signals in actuation-blocks into
    %processing-blocks. The program assumes that all the signals to be
    %reshaped are organized in the same manner. Inputs need to be in the
    %top level, i.e. nf (rather than nf.ablocks.smp). This program will
    %automatically work on .smp only (since .p makes no sense).
    
    %Get Organization Constants
    S = size(p_p);
    NPp = pp.sBS;
    NB = length(a_i); %number of raw processing blocks
    [NPa,NC,~] = size(sm_p_a); %size of actuation points in subblocks and number of channels
    RtS = S(3)/NB; %raw to sub-block conversion
    nNBps = floor(NPa*length(a_i{1})/NPp)*NB; %number of smoothed processing sub-blocks
    StS = nNBps/NB; %Smoothed to sub-block conversion
    
    %Initialize outputs
    output = zeros(NPp,NC,nNBps);
    
    %Replace blocks
    ablocknum = 1;
    pblocknum = 1;
    for n = 1:NB
        sig = p_p(:,:,n);
        
        for k = 1:length(a_i{n})
            sig(a_i{n}(k)+(0:NPa-1),:) = sm_p_a(:,:,ablocknum);
            ablocknum = ablocknum +1;
        end
        gb = rem(NPp*RtS-a_i{n}(1),NPp)+1;
        sig = sig(a_i{n}(1):end-gb,:);
        sig = reshape(sig,NPp,[],NC);
        output(:,:,pblocknum+(0:StS-1)) = permute(sig,[1 3 2]);
        pblocknum = pblocknum+StS;
    end
end

function SaveData(data,fname)
    %unpack structure
    nf = data.nf;
    ff = data.ff;
    phys = data.phys;
    src_dir = data.src_dir;
    filename = data.filename;
    trigger = data.trigger;
    test = data.test;

    save(fname,'nf','ff','phys','src_dir','filename','trigger','test');
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

function [output a_indices] = iSubBlocks(signal,trigger,pp)
    %Identifies forcing indices and reshapes signal into forcing-blocks,
    %rather than acquisition blocks
    
    %Identify actuation indices for each block
    a_indices = cell(1,pp.NB);
    t = cell(1,pp.NB);
    for n = 1:pp.NB
        a_indices{n} = IdentifyActuation(trigger(:,n)); %actuation indices
        t{n} = diff(a_indices{n});
    end
    T = round(mean(vertcat(t{:}))); %period of actuation, in indices 
    
    %Remove incomplete periods that occur at end of blocks
    for n = 1:pp.NB
        a_indices{n} = round(a_indices{n});
        a_indices{n} = a_indices{n}(a_indices{n}+T-1 < pp.BS);
    end
    
    %Reshape signal
    NSB = sum(cellfun(@(x) length(x),a_indices)); %total number of sub-blocks
    [~,M,~] = size(signal);
    output = zeros(T,M,NSB);
    counter = 1;
    for n = 1:pp.NB
        for nn = 1:length(a_indices{n})
            output(:,:,counter) = signal(a_indices{n}(nn):a_indices{n}(nn)+T-1,:,n);
            counter = counter + 1;
        end
    end

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
    actuation = actuation+double(trigger(actuation-1) > 0); %shifts indices by one if the found peak is on the voltage drop, rather than the voltage rise
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

function [sm_wave] = wavefilter(avg_wvfm,phys,pp)
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
    hww = 17;%15; %hard wavelet smoothing window
    sww = 9;%7; %soft wavelet smoothing window
    htw = 11;%15; %hard temporal smoothing window
    stw = 3;%7; %soft temporal smoothing window
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
