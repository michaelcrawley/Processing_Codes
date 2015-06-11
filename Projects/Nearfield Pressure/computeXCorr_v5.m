function [cormap lags maxcor maxlags X Y] = computeXCorr_v5(src_dir,flist,FFCh,flag)

    %this function computes the cross-correlation values between the
    %filtered, phase-averaged, near-field microphones to the specified
    %far-field microphone. The last input (flag) tells the program whether
    %to use the raw signal 'p' or filtered signal 'smp'. By default it uses
    %the filtered signal.
    
    if ~exist('flag','var'); flag = 'smp'; end
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 1:8;
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
    
    h = waitbar(0,'Calculating...');
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = eval(['data.ff.ablocks.',flag,'(:,pp.FFCh,:);']);
        nf = eval(['data.nf.ablocks.',flag,'(:,pp.NFCh,:);']);
        [sBS,~,NsBS] = size(nf);
                
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        for nn = 1:Lc
            cormap{n,nn} = zeros(sBS*2-1,1);
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(ff(:,1,nnn),nf(:,nn,nnn),'coeff'); 
               
               cormap{n,nn} = cormap{n,nn}+c/NsBS;
               lags{n,nn} = ilags/pp.FS;

            end
            [maxcor(n,nn),I] = max(cormap{n,nn},[],1);
            maxlags(n,nn) = lags{n,nn}(I);
        end
        waitbar(n/Lf,h,['Calculating...',num2str(round(n/Lf*100)),'% Complete']);
    end
    close(h);
    
    %Reshape and Sort spatial grid
    [X Y maxcor maxlags cormap lags] = ReshapeGrid(X,Y,maxcor,maxlags,cormap,lags);    
end