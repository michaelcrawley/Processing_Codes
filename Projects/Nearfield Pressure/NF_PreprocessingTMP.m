function NF_PreprocessingTMP(src_dir,flist,cal_dir,out_dir,cpus)
    %This is a test function, the purpose of which is to preprocess the 
    %data files. For the forced cases, it will filter out the
    %actuator self noise of the nearfield microphone signals before
    %phase-averaging (as opposed to the previous method where it was done
    %afterwards). 
    %Code version: 1.0
    %Output mat-file version: 1.1
    
    %Temporary Constants
    pp.BS = 81920; %block size of data
    pp.sBS = 8192; %sub-block size
    pp.NB = 10;
    pp.sNBmax = 100;
    pp.NFCh = 13:20; %nearfield microphone channels
    pp.lNFCh = 12:19; %nearfield channels after trigger has been removed
    pp.FFCh = 1:11; %farfield microphone channels    
    pp.tCh = 12; %channel number of trigger signal
    pp.NCh = length([pp.NFCh pp.FFCh])+1; %total number of channels recorded
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.D = 0.0254;
    pp.AT = 26.1; %atmospheric temperature, in C
    pp.sp = 1;
    pp.nm = 8;
    pp.R = [2.45 2.57 2.69 2.90 3.19 3.46 3.68 3.28 3.15 2.63 2.56]; %radial distance of farfield microphones, in meters
    
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
        tl{n}.phys = NF_ReadFileName(fname,pp);
        
        %grab processing info
        tl{n}.src_dir = src_dir;
        tl{n}.filename = filename;
        
        %Read in data 
        [nfp,ffp,trigger] = ReadData(src_dir,flist{n},pp,tl{n}.phys,cal_dir); 
                        
        %reshape into processing sub-blocks
        tl{n}.nf.pblocks.p = permute(nfp, [1 3 2]);
        tl{n}.ff.pblocks.p = permute(ffp, [1 3 2]);
        tl{n}.nf.pblocks.p = reshape(tl{n}.nf.pblocks.p,pp.sBS,[],length(pp.NFCh));
        tl{n}.ff.pblocks.p = reshape(tl{n}.ff.pblocks.p,pp.sBS,[],length(pp.FFCh));
        tl{n}.nf.pblocks.p = permute(tl{n}.nf.pblocks.p, [ 1 3 2]);        
        tl{n}.ff.pblocks.p = permute(tl{n}.ff.pblocks.p, [ 1 3 2]);
        
        if tl{n}.phys.Stdf ~= 0  %forced cases
            %Reshape data into forcing sub-blocks
            [pressure a_i]= iSubBlocks([ffp nfp],trigger,pp);           

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

            %need to figure out how to correctly break into processing
            %sub-blocks!
            slocator = 1;
            locator = 1;
            sm_pressure = permute(sm_pressure,[1 3 2]);
            smp = zeros(pp.sBS,pp.sNBmax,pp.NCh-1);
            for nn = 1:pp.NB
                nsblocks = length(a_i{nn});
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
        end 
        SaveData(tl{n},[out_dir filesep fname '_TMPproc.mat']);
        tl{n} = [];
    end
    matlabpool('close');
end

function SaveData(data,fname)
    %unpack structure
    nf = data.nf;
    ff = data.ff;
    phys = data.phys;
    src_dir = data.src_dir;
    filename = data.filename;

    save(fname,'nf','ff','phys','src_dir','filename');
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
            calfile = getfiles(['CAL*Ch' num2str(n) '*mVPa' num2str(phys.gain(n)) '*.mat'],cal_dir);
            tmp = load([cal_dir filesep calfile{1}]);
            cal(n) = tmp.PaV;
        end
        CAL = spdiags(cal,0,M,M);
        
        %Calibrate signals
        c = s*CAL; %convert pressure from volts to Pa 
        
        %Reshape back into Points x Channels x Blocks
        c = permute(reshape(c,N,O,M),[1 3 2]);
end

function [phys] = NF_ReadFileName(filename,pp)
%Parses specified filename to determine jet operating conditions (TTR, exit
%velocity, acoustic Mach number, forcing frequency, etc).  Fields are added
%to the 'phys' structure.
    
    md = getMetaFromStr(filename,'NearFieldParams'); %pull metadata from filename
    
    %pull necessary data
    M = md.M.value;
    T = md.T.value;
    Stdf = md.S.value;
    x = (md.x.value:pp.sp:md.x.value+(pp.nm-1)*pp.sp); %axial distances of microphones
    y = md.r.value+(x-x(1))*tand(md.a.value); % y/D
    r = [pp.R sqrt(x.^2+y.^2)*0.0254];%radial distance of nearfield microphones, converted to meters
    gain = [md.mVf.value*ones(length(pp.FFCh),1); 0]; %get farfield microphone gains, add in zero for trigger.  Assumes farfield mics were at same gain, and trigger is between farfield and nearfield mics
    if isfield(md,'mV')
        gain = [gain; md.mV.value*ones(pp.nm,1)];
    elseif isfield(md,'mVu') && isfield(md,'mVd');
        gain = [gain; md.mVu.value*ones(pp.nm/2,1); md.mVd.value*ones(pp.nm/2,1)]; %assumes gain settings were split in half over the NF array
    end
    
    
    %Calc jet/ambient properties
    To = T + 273.15; %convert to kelvin
    Te = To/(1+1/5*M^2); %jet exit temperature, in kelvin
    TTR = To/(pp.AT+273.15); %jet temperature ratio
    Ue = M*sqrt(1.4*287.05*Te); %nozzle exit velocity (m/s)
    Ma = M*sqrt(Te/(pp.AT+273.15)); %acoustic jet Mach number
    a = sqrt(1.4*287.05*(pp.AT+273.15)); %ambient speed of sound (m/s)
    c = Ue/M; %speed of sound in jet (m/s)
    Uc = Ue*a/(a+c); %estimated convective velocity
    t_c = r/Uc; %convective time for structures
    i_c = round(t_c*pp.FS); %convective index for structures
    t_a = r/a; %acoustic time for actuator pulse
    i_a = round(t_a*pp.FS); %acoustic index for actuator pulse
    
    phys = struct('M',M,'To',To,'x',x,'y',y,'r',r,'TTR',TTR,'Ue',Ue,'Ma',Ma,'a',a,'c',c,'Uc',Uc,'t_c',t_c,'i_c',i_c,'t_a',t_a,'i_a',i_a,'gain',gain,'Stdf',Stdf);
end

function [output a_indices] = iSubBlocks(signal,trigger,pp)
    %Identifies forcing indices and reshapes signal into forcing-blocks,
    %rather than acquisition blocks
    
    %Identify actuation indices for each block
    a_indices = cell(1,pp.NB);
    t = zeros(1,pp.NB);
    for n = 1:pp.NB
        a_indices{n} = round(IdentifyActuation(trigger(:,n))); %actuation indices
        t(n) = mean(diff(a_indices{n}));
    end
    T = round(mean(t)); %period of actuation, in indices
    
    %Remove incomplete periods that occur at end of blocks
    for n = 1:10
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
    hww = 19;%15; %hard wavelet smoothing window
    sww = 11;%7; %soft wavelet smoothing window
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
        
        %find actuator self-noise window
        tw = [center-wt center+wt]; %temporal index of window 
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
