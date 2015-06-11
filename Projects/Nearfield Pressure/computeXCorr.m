function [map maxcor lags X Y] = computeXCorr(pf,FFCh)
    %this function computes the cross-correlation values between the
    %phase-averaged, filtered, near-field microphones to the specified
    %far-field microphone.

    %Constants
    NFCh = 12:19; %nearfield microphone channel
%     FFCh = 7; %farfield microphone channels  
    FS = 2e5; %data acquisition frequency, in Hz

    %Initialize variables
    Lf = length(pf.sm_wvfm);
    Lc = length(NFCh);
    map.cor = cell(Lf,1);
    map.lags = cell(Lf,1);
    maxcor = zeros(Lf,Lc);
    lags = maxcor;
    X = maxcor;
    Y = maxcor;
    
    for n = 1:Lf
        X(n,:) = pf.phys{n}.x;
        Y(n,:) = pf.phys{n}.y;
        for nn = NFCh
           %compute nearfield to farfield correlations
           [c,ilags] = xcorr(pf.sm_wvfm{n}(:,FFCh),pf.sm_wvfm{n}(:,nn),'coeff'); 
           i = nn-NFCh(1)+1;
           map.cor{n}(:,i) = c;
           map.lags{n} = ilags/FS;
           [maxcor(n,i),I] = max(c);
           lags(n,i) = map.lags{n}(I); 
        end        
    end
    
    %Reshape and Sort spatial grid
    [X Y maxcor lags] = ReshapeGrid(X,Y,maxcor,lags);
end