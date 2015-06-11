function [cormap lags maxcor maxlags X Y] = computeNFXcorr(src_dir,flist,flag)

    if ~exist('flag','var'); flag = 'smp'; end
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 1:8; 

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    cormap = cell(Lf,Lc-1);
    lags = cormap;
    maxcor = zeros(Lf,Lc-1);
    maxlags = maxcor;
    X = zeros(Lf,Lc);
    Y = X;
    
    h = waitbar(0,'Calculating...');
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        nf = eval(['data.nf.pblocks.',flag,'(:,pp.NFCh,:);']);
        [sBS,~,NsBS] = size(nf);
                
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        for nn = 1:Lc-1
            cormap{n,nn} = zeros(sBS*2-1,1);
            for nnn = 1:NsBS
               %compute nearfield to nearfield correlations
               [c,ilags] = xcorr(nf(:,nn+1,nnn),nf(:,nn,nnn),'coeff');                
               cormap{n,nn} = cormap{n,nn}+c/NsBS;
            end
            lags{n,nn} = ilags/pp.FS;
            [maxcor(n,nn),I] = max(cormap{n,nn},[],1);
            maxlags(n,nn) = lags{n,nn}(I);
        end
        waitbar(n/Lf,h,['Calculating...',num2str(round(n/Lf*100)),'% Complete']);
    end
    close(h);    
end

