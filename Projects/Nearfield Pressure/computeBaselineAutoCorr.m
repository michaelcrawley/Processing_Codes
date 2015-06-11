function [nfcormap ffcormap lags X Y] = computeBaselineAutoCorr(src_dir,flist)

    %Temporary constants
    pp.BS = 81920;
    pp.sBS = 8192;
    pp.NB = 10;
    pp.NCh = 20;
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.FFCh = 1:11;
    pp.NFCh = 13:20;
    pp.D = 0.0254;
    pp.AT = 26.1; %atmospheric temperature, in C
    pp.sp = 1;
    pp.nm = 8;
    pp.R = 2.57; %radial distance of farfield microphones, in meters   

    %Initialize variables
    Lf = length(flist);
    Lc = pp.NCh-1; %remove trigger
    cormap = zeros(Lf,Lc,2*pp.sBS-1);
    X = zeros(Lf,length(pp.NFCh));
    Y = X;
    
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
        data = reshape(data,pp.sBS,[],pp.NCh);
        [~,NsBS,~] = size(data);
        
        channel = 0;
        for nn = [pp.FFCh pp.NFCh]
            channel = channel+1;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(data(:,nnn,nn),'coeff'); 
               
               cormap(n,channel,:) = squeeze(cormap(n,channel,:))+c/NsBS;
               lags = ilags/pp.FS;
            end
        end        
    end
    
    ffcormap = cormap(:,pp.FFCh,:);
    nfcormap = cormap(:,pp.NFCh-1,:);
    
    %Sort spatial grid
    %sort first so X values are ascending order
    [~,I] = sort(X(:,1));
    X = X(I,:);
    Y = Y(I,:);
    nfcormap = nfcormap(I,:,:);

    %Reshape matrix
    [~,I] = unique(X(:,1),'first');
    X = [X(I(1):I(2)-1,:) X(I(2):end,:)];
    Y = [Y(I(1):I(2)-1,:) Y(I(2):end,:)];
    nfcormap = [nfcormap(I(1):I(2)-1,:,:) nfcormap(I(2):end,:,:)];

    %Sort in X direction
    [~,I] = sort(X(1,:));
    X = X(:,I);
    Y = Y(:,I); 
    nfcormap = nfcormap(:,I,:);

    %Sort in Y direction
    [~,I] = sort(Y(:,1));
    X(:,1:2:end) = X(I,1:2:end);
    Y(:,1:2:end) = Y(I,1:2:end);
    nfcormap(:,1:2:end,:) = nfcormap(I,1:2:end,:);
    [~,I] = sort(Y(:,2));
    X(:,2:2:end) = X(I,2:2:end);
    Y(:,2:2:end) = Y(I,2:2:end);
    nfcormap(:,2:2:end,:) = nfcormap(I,2:2:end,:);
end