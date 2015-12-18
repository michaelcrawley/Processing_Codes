function [X Y lags TTcormap TPcormap TNcormap PPcormap NNcormap PNcormap NPcormap] = computeXCorr_v7(src_dir,flist,FFCh,flag)

    %this function computes the cross-correlation values between the
    %filtered, phase-averaged, near-field microphones to the specified
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
    TTcormap = cell(Lf,Lc);
    TPcormap = TTcormap;
    TNcormap = TTcormap;
    PPcormap = TTcormap;
    NNcormap = TTcormap;
    PNcormap = TTcormap;
    NPcormap = TTcormap;
    lags = TTcormap;
    
    for n = 1:Lf
        %Read file
        data = load([src_dir filesep flist{n}]);
        ff = mean(eval(['data.ff.ablocks.',flag,'(:,pp.FFCh,:);']),3);
        nf = mean(eval(['data.nf.ablocks.',flag,'(:,pp.NFCh,:);']),3);
        
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
            TTcormap{n,nn} = xcorr(ff,nf(:,nn),'coeff');                    
        end    
        
        %Positive Farfield - Total Nearfield Waveform
        for nn = 1:Lc
            %compute nearfield to farfield correlations
            TPcormap{n,nn} = xcorr(Pff,nf(:,nn),'coeff');                
        end  
        
        %Positive Farfield - Positive Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            PPcormap{n,nn} = xcorr(Pff,Pnf,'coeff');                
        end  
        
        %Negative Farfield - Total Nearfield Waveform
        for nn = 1:Lc
            %compute nearfield to farfield correlations
            TNcormap{n,nn} = xcorr(Nff,nf(:,nn),'coeff');               
        end 

        %Negative Farfield - Negative Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            NNcormap{n,nn} = xcorr(Nff,Nnf,'coeff');                
        end 
        
        %Positive Farfield - Negative Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
               %compute nearfield to farfield correlations
               PNcormap{n,nn} = xcorr(Nff,Pnf,'coeff');                
        end 
        
        %Negative Farfield - Positive Nearfield Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            %compute nearfield to farfield correlations
            [NPcormap{n,nn},ilags] = xcorr(Pff,Nnf,'coeff');                
            lags{n,nn} = ilags/pp.FS;
        end 
    end
    
    %Reshape and Sort spatial grid
    [X,Y,TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap,lags] = ReshapeGrid(X,Y,TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap,lags);
    [TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap,lags] = EqualizeCells(TTcormap,TPcormap,TNcormap,PPcormap,NNcormap,PNcormap,NPcormap,lags);
end