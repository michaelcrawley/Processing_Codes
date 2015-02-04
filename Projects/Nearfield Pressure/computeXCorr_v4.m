function [cormap lags maxcor maxlags X Y] = computeXCorr_v4(src_dir,flist,FFCh,flag)

    %this function computes the cross-correlation values between the
    %filtered, non-phase-averaged, near-field microphones to the specified
    %far-field microphone. The last input (flag) tells the program whether
    %to use the raw signal 'p' or filtered signal 'smp'. By default it uses
    %the filtered signal.
    
    if ~exist('flag','var'); flag = 'smp'; end
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 1:16;
    pp.FFCh = FFCh;
    pp.sBS = 8192;

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    cormap = zeros(Lf,Lc,2*pp.sBS-1);
    maxcor = zeros(Lf,Lc);
    maxlags = maxcor;
    X = maxcor;
    Y = maxcor;
    
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = eval(['data.ff.pblocks.',flag,'(:,pp.FFCh,:);']);
        nf = eval(['data.nf.pblocks.',flag,'(:,pp.NFCh,:);']);
        [~,~,NsBS] = size(nf);
                
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(ff(:,1,nnn),nf(:,nn,nnn),'coeff'); 
               
               cormap(n,nn,:) = squeeze(cormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [maxcor(n,nn),I] = max(cormap(n,nn,:),[],3);
            maxlags(n,nn) = lags(I);
        end        
    end
    
    %Reshape and Sort spatial grid
    [X Y maxcor maxlags cormap] = ReshapeGrid(X,Y,maxcor,maxlags,cormap);    
end