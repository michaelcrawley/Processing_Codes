function [cormap lags maxcor maxlags X Y] = computeBaselineXCorr(src_dir,flist,FFCh)

    %Temporary constants
    pp.BS = 81920;
    pp.sBS = 8192;
    pp.NB = 10;
    pp.NCh = 20;
    pp.NFCh = 13:20; %nearfield microphone channel
    pp.FFCh = FFCh; %farfield microphone channels for correlation
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.D = 0.0254;
    pp.AT = 26.1; %atmospheric temperature, in C
    pp.sp = 1;
    pp.nm = 8;
    pp.R = 2.57; %radial distance of farfield microphones, in meters   

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    cormap = zeros(Lf,Lc,2*pp.sBS-1);
    maxcor = zeros(Lf,Lc);
    maxlags = maxcor;
    X = maxcor;
    Y = maxcor;
    
    for n = 1:Lf
        %determine x and y locations for microphone signals
        [~,fname] = fileparts(flist{n});
        phys = NF_ReadFileName(fname,pp);
        X(n,:) = phys.x;
        Y(n,:) = phys.y;
       
        %Read file
        fid = fopen([src_dir filesep flist{n}],'r');
        raw = -fread(fid,'float32');fclose(fid);
        data = reshape(raw,pp.BS,pp.NCh,[]); %reshape data into Points x Channels x Blocks
        data = permute(data,[1 3 2]); %reshape data into Points x Blocks x Channels
        ff = reshape(data(:,:,pp.FFCh),pp.sBS,[],length(pp.FFCh)); 
        nf = reshape(data(:,:,pp.NFCh),pp.sBS,[],length(pp.NFCh));
        [~,NsBS] = size(ff);
        
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(ff(:,nnn),nf(:,nnn,nn),'coeff'); 
               
               cormap(n,nn,:) = squeeze(cormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
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