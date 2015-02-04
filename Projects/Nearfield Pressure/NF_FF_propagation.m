function [prop_sig,r,R] = NF_FF_propagation(sig,x,xs,xf,units,a,FS)
%Propagates near-field signal to far-field position, assuming spherical
%(acoustic) spreading. 
%Inputs:
%           sig:        signal to be propagated (NP x NC)
%           x:          location of near-field microphone (NC x 1)
%           xs:         assumed source location [x,y]
%           xf:         location of far-field microphone
%           units:      'pa', 'pa^2', or 'dB', depending on the format of
%                       the signal
%           a:          ambient speed of sound, only used for time-domain
%                       propagation of signal (units must be 'pa')
%           FS:         sampling frequency (Hz), only used for time-domain
%                       propagation
    
    %Calculate distances
    r = sqrt((x(:,1)-xs(1)).^2 + (x(:,2)-xs(2)).^2); %distance from source to near-field
    R = sqrt(sum((xf-xs).^2)); %distance from source to far-field
    ratio = r(:)'/R;
    
    [NP,NC] = size(sig);
    switch lower(units)
        case 'pa'
            offset = ratio;
            prop_sig = sig.*repmat(offset,NP,1);
            
            %calculate propagation times
            i_n = round((r/a)*FS); %indices from source to near-field
            i_f = round((R/a)*FS); %indices from source to far-field
            for n = 1:NC
                prop_sig(:,n) = circshift(prop_sig(:,n),[(i_f-i_n(n)),0]);
            end
        case 'pa^2'
            offset = ratio.^2;
            prop_sig = sig.*repmat(offset,NP,1);
        case 'db'
            offset = 20*log10(ratio);
            prop_sig = sig + repmat(offset,NP,1);
        otherwise
            error('Wrong units supplied');
    end

end