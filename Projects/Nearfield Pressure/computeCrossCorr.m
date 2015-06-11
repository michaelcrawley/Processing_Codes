function computeCrossCorr(src_dir,flist,s_ref,s_sig,savename,pp,flags)

    %Note that s_sig can only contain near-field microphones
    
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
    
    %Parse Signal & Reference mic commands
    Lcp = length(s_ref); %Number of reference-signal pairs
    split_sig = regexp(s_sig,' ','split');
    split_ref = regexp(s_ref,' ','split');
    exp = '^(?<tag>[nf]f)?(?<txt>[0-9\:]+)';
    for n = 1:Lcp
        cmds_sig{n} = regexpi(split_sig{n},exp,'names');
        cmds_ref{n} = regexpi(split_ref{n},exp,'names');
        Lcs(n) = 0;
        Lcr(n) = 0;
        sig_chn{n} = [];
        for q = 1:length(cmds_sig{n})
            cmds_sig{n}{q}.num = eval(cmds_sig{n}{q}.txt);
            sig_chn{n} = [sig_chn{n} cmds_sig{n}{q}.num]; 
            Lcs(n)  = Lcs(n)+length(cmds_sig{n}{q}.num);
        end
        for q = 1:length(cmds_ref{n})
            cmds_ref{n}{q}.num = eval(cmds_ref{n}{q}.txt);
            Lcr(n)  = Lcr(n)+length(cmds_ref{n}{q}.num);  %#ok<*AGROW>
        end
    end
    
    %Initialize variables
    Lf = length(flist);
    cormap = cell(Lcp,1);
    X = cormap;
    Y = X;
    for n = 1:Lcp
        X{n} = zeros(Lf,Lcs(n));
        Y{n} = X{n};
        if cellflag
            lags{n} = cell(Lf,1);
            cormap{n} = cell(Lf,Lcs(n),Lcr(n));
        else
            cormap{n} = zeros(Lf,Lcs(n),2*pp.sBS-1,Lcr(n));
        end
    end

    for n = 1:Lf %All Files
        %Read file
        data = load([src_dir filesep flist{n}]);
        
        for m = 1:Lcp %All signal-reference pairs
            sig = [];
            for q = 1:length(cmds_sig{m})
                sig = [sig eval(['data.',cmds_sig{m}{q}.tag,'.',procflag,'(:,',cmds_sig{m}{q}.txt,',:);'])];
            end
            ref = [];
            for q = 1:length(cmds_ref{m})
                ref = [ref eval(['data.',cmds_ref{m}{q}.tag,'.',procflag,'(:,',cmds_ref{m}{q}.txt,',:);'])];
            end
            [BS,Nsig,NsBS] = size(sig); %Get data block size and number of (sub) blocks

            %determine x and y locations for microphone signals
            X{m}(n,:) = data.phys.x(sig_chn{m});
            Y{m}(n,:) = data.phys.y(sig_chn{m});

            %Replicate signal, if desired (for phase-averaged waveforms)
            if repflag
                Nrep = ceil(pp.sBS/BS);
                BS = pp.sBS;
                if avgflag  %phase-average waveforms
                    sig = mean(sig,3);
                    ref = mean(ref,3);
                    sig = repmat(sig,Nrep,1);
                    sig = sig(1:BS,:,:);
                    ref = repmat(ref,Nrep,1);
                    ref = ref(1:BS,:,:); 
                    NsBS = 1;
                else %create swap variables if the reps will not be done to phaavg
                    sig_swap = sig;
                    ref_swap = ref;
                end            
            end
                        
            for r = 1:Lcr(m)
                %Compute Cross-Correlations
                for q = 1:NsBS
                    if repflag && ~avgflag %replicate each block - needed for memory issues
                        sig = sig_swap(:,:,q);
                        ref = ref_swap(:,:,q);
                        sig = repmat(sig,Nrep,1);
                        sig = sig(1:BS,:,:);
                        ref = repmat(ref,Nrep,1);
                        ref = ref(1:BS,:,:); 
                    end
                    for k = 1:Nsig
                        [c,ilags] = xcorr(ref(:,r,q),sig(:,k,q),'coeff');
                        if cellflag
                            cormap{m}{n,k,r} = c;
                        else
                            cormap{m}(n,k,:,r) = squeeze(cormap{m}(n,k,:,r))+c/NsBS;
                        end
                    end
                end
            end
            if cellflag
                lags{m}{n} = ilags/pp.FS;
            end
        end
    end
    if ~cellflag
        lags = ilags/pp.FS;
    end
        
    %Reshape and Sort spatial grid
    for m = 1:Lcp
        if cellflag
            [~,~,lags{m}] = ReshapeGrid(X{m}(:,1),Y{m}(:,1),lags{m});
        end
        [X{m} Y{m} cormap{m}] = ReshapeGrid(X{m},Y{m},cormap{m}); 
    end
    
    %Save Outputs
    save(savename,'X','Y','cormap','lags','pp','procflag','repflag','s_ref','s_sig');
end