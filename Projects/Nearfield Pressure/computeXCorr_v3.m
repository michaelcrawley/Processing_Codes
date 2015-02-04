function [cormap lags maxcor maxlags X Y] = computeXCorr_v3(src_dir,flist)

    %this function computes the cross-correlation values between the
    %filtered, non-phase-averaged, near-field microphones to the specified
    %far-field microphone.
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 2:9;
    pp.FFCh = 1;
    pp.sBS = 8192;

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    cormap = zeros(Lf,Lc,pp.sBS*2-1);
    maxcor = zeros(Lf,Lc);
    maxlags = maxcor;
    X = maxcor;
    Y = maxcor;
    
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        
        %Reshape file
        pressure = reshape(permute(data.sm_pressure,[1 3 2]),[],9);
        ff = pressure(:,pp.FFCh);
        nf = pressure(:,pp.NFCh);
        gb = rem(length(ff),pp.sBS);
        ff = reshape(ff(1:end-gb),pp.sBS,[]);
        nf = permute(reshape(nf(1:end-gb,:),pp.sBS,[],length(pp.NFCh)),[1 3 2]);
        [~,NsBS] = size(ff);
                
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(ff(:,nnn),nf(:,nn,nnn),'coeff');                
               cormap(n,nn,:) = squeeze(cormap(n,nn,:))+c/NsBS;
            end
            lags = ilags/pp.FS;
            [maxcor(n,nn),I] = max(cormap(n,nn,:),[],3);
            maxlags(n,nn) = lags(I);
        end        
    end
    
    %Sort spatial grid
    %sort first so X values are ascending order
    [~,I] = sort(X(:,1));
    X = X(I,:);
    Y = Y(I,:);
    maxcor = maxcor(I,:);
    maxlags = maxlags(I,:);
    cormap = cormap(I,:,:);

    %Reshape matrix
    [~,I] = unique(X(:,1),'first');
    X = [X(I(1):I(2)-1,:) X(I(2):end,:)];
    Y = [Y(I(1):I(2)-1,:) Y(I(2):end,:)];
    maxcor = [maxcor(I(1):I(2)-1,:) maxcor(I(2):end,:)];
    maxlags = [maxlags(I(1):I(2)-1,:) maxlags(I(2):end,:)];
    cormap = [cormap(I(1):I(2)-1,:,:) cormap(I(2):end,:,:)];

    %Sort in X direction
    [~,I] = sort(X(1,:));
    X = X(:,I);
    Y = Y(:,I); 
    maxcor = maxcor(:,I);
    maxlags = maxlags(:,I);
    cormap = cormap(:,I,:);

    %Sort in Y direction
    [~,I] = sort(Y(:,1));
    X(:,1:2:end) = X(I,1:2:end);
    Y(:,1:2:end) = Y(I,1:2:end);
    maxcor(:,1:2:end) = maxcor(I,1:2:end);
    maxlags(:,1:2:end) = maxlags(I,1:2:end);
    cormap(:,1:2:end,:) = cormap(I,1:2:end,:);
    [~,I] = sort(Y(:,2));
    X(:,2:2:end) = X(I,2:2:end);
    Y(:,2:2:end) = Y(I,2:2:end);
    maxcor(:,2:2:end) = maxcor(I,2:2:end);
    maxlags(:,2:2:end) = maxlags(I,2:2:end);
    cormap(:,2:2:end,:) = cormap(I,2:2:end,:);
    
end