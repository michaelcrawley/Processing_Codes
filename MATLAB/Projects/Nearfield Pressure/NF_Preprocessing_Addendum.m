function NF_Preprocessing_Addendum(proc_dir,src_dir,flist,params)

    %Constants
    load(params,'pp');
    pp.NCh = length([pp.NFCh pp.FFCh])+1;
    
    for n = 1:length(flist)
        data = load([proc_dir filesep flist{n}]);
        trigger.sig = ReadTrigger(src_dir,data.filename,pp);
        
        if data.phys.Stdf ~= 0  %forced cases
            trigger.a_i = iSubBlocks(trigger.sig,pp);
        end

        save([proc_dir filesep flist{n}],'-append','trigger','src_dir');
    end


end

function [trigger] = ReadTrigger(src_dir,fname,pp)
    %Read file
    fid = fopen([src_dir filesep fname],'r');
    raw = fread(fid,'float32');fclose(fid);
    tmp = reshape(raw,pp.BS,pp.NCh,[]); %reshape data into Points x Channels x Blocks

    %Extract components
    trigger = squeeze(tmp(:,pp.tCh,:)); %extract trigger from rest of voltage traces        
end

function [a_indices] = iSubBlocks(trigger,pp)
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