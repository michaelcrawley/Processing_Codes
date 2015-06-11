function [X Y lags TTcormap TPcormap TNcormap PPcormap NNcormap PNcormap NPcormap] = computeXCorr_v8(src_dir,flist,FFCh,flag)

    %this function computes the cross-correlation values between the
    %filtered, phase-averaged, replicated, near-field microphones to the specified
    %far-field microphone. Three correlations are actually computed: the
    %total far-field waveform, the positive excursions in the far-field
    %waveform, and the negative excursions in the far-field waveform.
    
    if ~exist('flag','var'); flag = 'smp'; end
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 1:8;
    pp.FFCh = FFCh;
    pp.BS = 8192;

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    X = zeros(Lf,Lc);
    Y = X;
    TTcormap = zeros(Lf,Lc,2*pp.BS-1);
    TPcormap = TTcormap;
    TNcormap = TTcormap;
    PPcormap = TTcormap;
    NNcormap = TTcormap;
    PNcormap = TTcormap;
    NPcormap = TTcormap;
    
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = mean(eval(['data.ff.ablocks.',flag,'(:,pp.FFCh,:);']),3);
        nf = mean(eval(['data.nf.ablocks.',flag,'(:,pp.NFCh,:);']),3);
        
        %replicate waveform
        Nrep = ceil(pp.BS/length(ff));
        ff = repmat(ff,Nrep,1);
        ff = ff(1:pp.BS,:);
        nf = repmat(nf,Nrep,1);
        nf = nf(1:pp.BS,:);
        
        meanff = mean(ff(:));
        Pff = zeros(size(ff));
        Pff(ff >= meanff) = ff(ff >= meanff);
        Nff = zeros(size(ff));
        Nff(ff <= meanff) = ff(ff <= meanff);
                
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        %Total Waveform
        for nn = 1:Lc
            %compute nearfield to farfield correlations
            c = xcorr(ff,nf(:,nn),'coeff');
            TTcormap(n,nn,:) = permute(c,[2 3 1]);
        end    
        
        %Positive Farfield - Total Nearfield Waveform
        for nn = 1:Lc
            %compute nearfield to farfield correlations
            c = xcorr(Pff,nf(:,nn),'coeff');
            TPcormap(n,nn,:) = permute(c,[2 3 1]);
        end  
        
        %Positive Farfield - Positive Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            c = xcorr(Pff,Pnf,'coeff'); 
            PPcormap(n,nn,:) = permute(c,[2 3 1]);
        end  
        
        %Negative Farfield - Total Nearfield Waveform
        for nn = 1:Lc
            %compute nearfield to farfield correlations
            c = xcorr(Nff,nf(:,nn),'coeff');
            TNcormap(n,nn,:) = permute(c,[2 3 1]);
        end 

        %Negative Farfield - Negative Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            c = xcorr(Nff,Nnf,'coeff'); 
            NNcormap(n,nn,:) = permute(c,[2 3 1]);
        end 
        
        %Positive Farfield - Negative Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            c = xcorr(Nff,Pnf,'coeff');   
            PNcormap(n,nn,:) = permute(c,[2 3 1]);
        end 
        
        %Negative Farfield - Positive Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            [c,ilags] = xcorr(Pff,Nnf,'coeff');                
            NPcormap(n,nn,:) = permute(c,[2 3 1]);
        end 
    end
    lags = ilags/pp.FS;
    
    %Reshape and Sort spatial grid
    [X,Y,TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap] = ReshapeGrid(X,Y,TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap);
end