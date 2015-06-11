function [cohmap f X Y phase nfpsd ffpsd] = computeCOH_v2(src_dir,flist)
%This function is designed to work on the filtered (non-averaged)
%microphone data. 


    %Temporary constants
    pp.BS = 81920;
    pp.sBS = 8192;
    pp.NB = 10;
    pp.NCh = 20;
    pp.NFCh = 2:9; %nearfield microphone channel
    pp.FFCh = 1; %farfield microphone channels  
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.D = 0.0254;
    pp.AT = 26.1; %atmospheric temperature, in C
    pp.sp = 1;
    pp.nm = 8;
    pp.R = 2.57; %radial distance of farfield microphones, in meters   
    
    %initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    nfrq = ceil((pp.sBS+1)/2);
    X = zeros(Lf,Lc);
    Y = X;
    cohmap = zeros(Lf,Lc,nfrq);
    phase = cohmap;
    nfpsd = cohmap;
    ffpsd = cohmap;
    
    for n = 1:Lf    
        %determine x and y locations for microphone signals
        [~,fname] = fileparts(flist{n});
        phys = NF_ReadFileName(fname,pp);
        X(n,:) = phys.x;
        Y(n,:) = phys.y;
       
        %Read file
        tmp = load([src_dir filesep flist{n}]);
        
        %Reshape file
        data = permute(tmp.sm_pressure,[1 3 2]);
        ff = data(:,:,pp.FFCh);        
        nf = data(:,:,pp.NFCh);
        ff = reshape(ff,[],length(pp.FFCh));
        nf = reshape(nf,[],length(pp.NFCh));
        
        %Cut file
        remainder = rem(length(ff),pp.sBS);
        ff = ff(1:end-remainder);
        nf = nf(1:end-remainder,:);
       
        %compute coherence
        for nn = 1:Lc
            [cohmap(n,nn,:), f, phase(n,nn,:), ~, nfpsd(n,nn,:), ffpsd(n,nn,:)] = coherence(nf(:,nn),ff,pp.FS,pp.sBS,@hann);          
        end
    end
    
    %sort first so X values are ascending order
    [~,I] = sort(X(:,1));
    X = X(I,:);
    Y = Y(I,:);
    cohmap = cohmap(I,:,:);
    phase = phase(I,:,:);
    nfpsd = nfpsd(I,:,:);
    ffpsd = ffpsd(I,:,:);

    %Reshape matrix
    [~,I] = unique(X(:,1),'first');
    X = [X(I(1):I(2)-1,:) X(I(2):end,:)];
    Y = [Y(I(1):I(2)-1,:) Y(I(2):end,:)];
    cohmap = [cohmap(I(1):I(2)-1,:,:) cohmap(I(2):end,:,:)];
    phase = [phase(I(1):I(2)-1,:,:) phase(I(2):end,:,:)];
    nfpsd = [nfpsd(I(1):I(2)-1,:,:) nfpsd(I(2):end,:,:)];
    ffpsd = [ffpsd(I(1):I(2)-1,:,:) ffpsd(I(2):end,:,:)];

    %Sort in X direction
    [~,I] = sort(X(1,:));
    X = X(:,I);
    Y = Y(:,I); 
    cohmap = cohmap(:,I,:);
    phase = phase(:,I,:);
    nfpsd = nfpsd(:,I,:);
    ffpsd = ffpsd(:,I,:);

    %Sort in Y direction
    [~,I] = sort(Y(:,1));
    X(:,1:2:end) = X(I,1:2:end);
    Y(:,1:2:end) = Y(I,1:2:end);
    cohmap(:,1:2:end,:) = cohmap(I,1:2:end,:);
    phase(:,1:2:end,:) = phase(I,1:2:end,:);
    nfpsd(:,1:2:end,:) = nfpsd(I,1:2:end,:);
    ffpsd(:,1:2:end,:) = ffpsd(I,1:2:end,:);    
    [~,I] = sort(Y(:,2));
    X(:,2:2:end) = X(I,2:2:end);
    Y(:,2:2:end) = Y(I,2:2:end);
    cohmap(:,2:2:end,:) = cohmap(I,2:2:end,:);
    phase(:,2:2:end,:) = phase(I,2:2:end,:);
    nfpsd(:,2:2:end,:) = nfpsd(I,2:2:end,:);
    ffpsd(:,2:2:end,:) = ffpsd(I,2:2:end,:);
end