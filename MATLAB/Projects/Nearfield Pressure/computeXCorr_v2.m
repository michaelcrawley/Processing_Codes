function [cormap lags maxcor maxlags X Y] = computeXCorr_v2(src_dir,flist,FFCh)

    %this function computes the cross-correlation values between the
    %filtered, non-phase-averaged, near-field microphones to the specified
    %far-field microphone.
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 2:9;
    pp.FFCh = FFCh;

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    cormap = cell(Lf,Lc);
    lags = cormap;
    maxcor = zeros(Lf,Lc);
    maxlags = maxcor;
    X = maxcor;
    Y = maxcor;
    
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = data.sm_pressure(:,pp.FFCh,:);
        nf = data.sm_pressure(:,pp.NFCh,:);
        [Ld,~,NsBS] = size(data.sm_pressure);
                
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        for nn = 1:Lc
            cormap{n,nn} = zeros(Ld*2-1,1);
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(ff(:,1,nnn),nf(:,nn,nnn),'coeff'); 
               
               cormap{n,nn} = cormap{n,nn}+c/NsBS;
               lags{n,nn} = ilags/pp.FS;

            end
            [maxcor(n,nn),I] = max(cormap{n,nn});
            maxlags(n,nn) = lags{n,nn}(I);
        end        
    end
    
    %Reshape and Sort spatial grid
    [X Y maxcor maxlags cormap] = ReshapeGrid(X,Y,maxcor,maxlags,cormap);    
end