function [X Y lags TTcormap TPcormap TNcormap PPcormap NNcormap PNcormap NPcormap] = computeXCorr_v6(src_dir,flist,FFCh,flag)

    %this function computes the cross-correlation values between the
    %filtered, non-phase-averaged, near-field microphones to the specified
    %far-field microphone. Three correlations are actually computed: the
    %total far-field waveform, the positive excursions in the far-field
    %waveform, and the negative excursions in the far-field waveform.
    
    if ~exist('flag','var'); flag = 'smp'; end
    
    %Temporary constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.NFCh = 1:8;
    pp.FFCh = FFCh;
    pp.sBS = 8192;

    %Initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    X = zeros(Lf,Lc);
    Y = X;
    TTcormap = zeros(Lf,Lc,2*pp.sBS-1);
    TPcormap = TTcormap;
    TNcormap = TTcormap;
    PPcormap = TTcormap;
    NNcormap = TTcormap;
    PNcormap = TTcormap;
    NPcormap = TTcormap;
    
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = eval(['data.ff.pblocks.',flag,'(:,pp.FFCh,:);']);
        nf = eval(['data.nf.pblocks.',flag,'(:,pp.NFCh,:);']);
        [~,~,NsBS] = size(nf);
        
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
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(ff(:,1,nnn),nf(:,nn,nnn),'coeff');                
               TTcormap(n,nn,:) = squeeze(TTcormap(n,nn,:))+c/NsBS;
            end
        end    
        
        %Positive Farfield - Total Nearfield Waveform
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Pff(:,1,nnn),nf(:,nn,nnn),'coeff');                
               TPcormap(n,nn,:) = squeeze(TPcormap(n,nn,:))+c/NsBS;
            end
        end  
        
        %Positive Farfield - Positive Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Pff(:,1,nnn),Pnf(:,nnn),'coeff');                
               PPcormap(n,nn,:) = squeeze(PPcormap(n,nn,:))+c/NsBS;
            end
        end  
        
        %Negative Farfield - Total Nearfield Waveform
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Nff(:,1,nnn),nf(:,nn,nnn),'coeff');               
               TNcormap(n,nn,:) = squeeze(TNcormap(n,nn,:))+c/NsBS;
            end
        end 

        %Negative Farfield - Negative Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Nff(:,1,nnn),Nnf(:,nnn),'coeff');                
               NNcormap(n,nn,:) = squeeze(NNcormap(n,nn,:))+c/NsBS;
            end
        end 
        
        %Positive Farfield - Negative Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Nff(:,1,nnn),Pnf(:,nnn),'coeff');                
               PNcormap(n,nn,:) = squeeze(PNcormap(n,nn,:))+c/NsBS;
            end
        end 
        
        %Negative Farfield - Positive Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Pff(:,1,nnn),Nnf(:,nnn),'coeff');                
               NPcormap(n,nn,:) = squeeze(NPcormap(n,nn,:))+c/NsBS;
            end
        end 
    end
    lags = ilags/pp.FS;
    
    %Reshape and Sort spatial grid
    [X,Y,TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap] = ReshapeGrid(X,Y,TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap);
end