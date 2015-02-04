function [peaks peak_locs] = IdentifyMaxPeaks(cor)

    npeaks = 2;
    [N,~] = size(cor);
    peaks = zeros(N,2);
    peak_locs = zeros(N,2);
    for n = 1:N
        [vals,loc] = findpeaks(cor(n,:),'minpeakheight',0);
        [vals,I] = sort(vals);
        loc = loc(I);
        
        %sort by amplitude
        vals = vals(end-npeaks+1:end);
        loc = loc(end-npeaks+1:end);
        
        %sort by time lag
        [peak_locs(n,:),I] = sort(loc);
        peaks(n,:) = vals(I);
    end

end