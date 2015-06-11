function [output a_indices] = iSubBlocks_EXP(signal,trigger,pp,offset)
    %Identifies forcing indices and reshapes signal into forcing-blocks,
    %rather than acquisition blocks. Nresample is the resample ratio (i.e.
    %for a physical sample rate of 200 kHz and Nresample = 10, the data is
    %resampled at 2 MHz). Nresample should be a power of 10. Norder is the
    %polynomial order for resampling.
    
    %Identify actuation indices for each block
    a_indices = cell(1,pp.NB);
    t = cell(1,pp.NB);
    for n = 1:pp.NB
        a_indices{n} = IdentifyActuation(trigger(:,n),pp,offset(n)); %actuation indices
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

function [a_i] = IdentifyActuation(trigger,pp,offset)
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
    
    %get rid of additional indices corresponding to channel firings which
    %are NOT channel 1
    a_i = a_i(offset:pp.NA:end);
end