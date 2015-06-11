function computeAutoCorr(src_dir,flist,savename,pp,flags)
    
    %Determine processing flags
    if ~exist('flags','var'), 
        procflag = 'pblocks.smp';
        repflag = false;
        avgflag = false;
        cellflag = false;
    else
        procflag = flags{1};
        repflag = any(strcmpi(flags,'rep')); %do the signals need to be repeated?
        avgflag = any(strcmpi(flags,'phavg')); %phase-average the waveforms before computing correlations?
        cellflag = strcmp('ablocks',procflag(1:7)) & ~repflag; %do the outputs need to be in cell array format?
    end
        
    %Initialize variables
    Lf = length(flist);
    Lcnf = length(pp.NFCh);
    Lcff = length(pp.FFCh);
    X = zeros(Lf,Lcnf);
    Y = X;
    if cellflag
        lags = cell(Lf,1);
        cormap.nf = cell(Lf,Lcnf);
        cormap.ff = cell(Lf,Lcff);
    else
        cormap.nf = zeros(Lf,Lcnf,2*pp.sBS-1);
        cormap.ff = zeros(Lf,Lcff,2*pp.sBS-1);
    end

    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = eval(['data.ff.',procflag,'(:,pp.FFCh,:);']);
        nf = eval(['data.nf.',procflag,'(:,pp.NFCh,:);']);
        [BS,~,NsBS] = size(nf); %Get data block size and number of (sub) blocks
        
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        %Replicate signal, if desired 
        if repflag 
            Nrep = ceil(pp.sBS/BS);
            BS = pp.sBS;
            if avgflag %phase-average waveforms
                nf = mean(nf,3);
                ff = mean(ff,3);
                nf = repmat(nf,Nrep,1);
                nf = nf(1:BS,:,:);
                ff = repmat(ff,Nrep,1);
                ff = ff(1:BS,:,:); 
                NsBS = 1;
            else %create swap variables if the reps will not be done to phaavg
                nf_swap = nf;
                ff_swap = ff;                
            end         
        end
        
        %Compute AutoCorrelations
        for q = 1:NsBS
            if repflag && ~avgflag %replicate each block - needed for memory issues
                nf = nf_swap(:,:,q);
                ff = ff_swap(:,:,q);                
                nf = repmat(nf,Nrep,1);
                nf = nf(1:BS,:,:);
                ff = repmat(ff,Nrep,1);
                ff = ff(1:BS,:,:);                
            end
            
            %Near-field Channels
            for k = pp.NFCh
                sig = squeeze(nf(:,k,q));
                sig = sig-mean(sig(:)); %normalize
                c = xcorr(sig,'coeff');
                if cellflag
                    cormap.nf{n,k} = c;
                else
                    cormap.nf(n,k,:) = squeeze(cormap.nf(n,k,:))+c/NsBS;
                end
            end
            
            %Far-field Channels
            for k = pp.FFCh
                sig = squeeze(ff(:,k,q));
                sig = sig-mean(sig(:)); %normalize
                [c,ilags] = xcorr(ff(:,k,q),'coeff');
                if cellflag
                    cormap.ff{n,k} = c;
                else
                    cormap.ff(n,k,:) = squeeze(cormap.ff(n,k,:))+c/NsBS;
                end
            end  
        end
        if cellflag
            lags{n} = ilags/pp.FS;
        end
    end
    if ~cellflag
        lags = ilags/pp.FS;
    end
        
    %Reshape and Sort spatial grid
    [~,~,cormap.ff] = ReshapeGrid(X(:,1:3),Y(:,1:3),cormap.ff);
    if cellflag
        [~,~,lags] = ReshapeGrid(X(:,1),Y(:,1),lags);
    end
    [X Y cormap.nf] = ReshapeGrid(X,Y,cormap.nf); %#ok<ASGLU>
    
    %Save Outputs
    save(savename,'X','Y','cormap','lags','pp','procflag','repflag','avgflag');
end