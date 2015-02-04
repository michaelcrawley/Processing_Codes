function [Tcormap lags Tmaxcor Tmaxlags X Y Pcormap Ncormap PPcormap NNcormap PNcormap] = computeXCorr_v6_old(src_dir,flist,FFCh,flag)

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
    Tcormap = zeros(Lf,Lc,2*pp.sBS-1);
    Tmaxcor = zeros(Lf,Lc);
    Tmaxlags = Tmaxcor;
    X = Tmaxcor;
    Y = Tmaxcor;
    Pcormap = zeros(Lf,Lc,2*pp.sBS-1);
    Pmaxcor = zeros(Lf,Lc);
    Pmaxlags = Tmaxcor;
    Ncormap = zeros(Lf,Lc,2*pp.sBS-1);
    Nmaxcor = zeros(Lf,Lc);
    Nmaxlags = Tmaxcor;
    PPcormap = Pcormap;
    PPmaxcor = Pmaxcor;
    PPmaxlags = Pmaxlags;
    NNcormap = Ncormap;
    NNmaxcor = Nmaxcor;
    NNmaxlags = Nmaxlags;
    PNcormap = Ncormap;
    PNmaxcor = Nmaxcor;
    PNmaxlags = Nmaxlags;
    
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
               
               Tcormap(n,nn,:) = squeeze(Tcormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [Tmaxcor(n,nn),I] = max(Tcormap(n,nn,:),[],3);
            Tmaxlags(n,nn) = lags(I);
        end    
        
        %Positive-Total Waveform
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Pff(:,1,nnn),nf(:,nn,nnn),'coeff'); 
               
               Pcormap(n,nn,:) = squeeze(Pcormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [Pmaxcor(n,nn),I] = max(Pcormap(n,nn,:),[],3);
            Pmaxlags(n,nn) = lags(I);
        end  
        
        %Positive-Positive Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Pff(:,1,nnn),Pnf(:,nnn),'coeff'); 
               
               PPcormap(n,nn,:) = squeeze(PPcormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [PPmaxcor(n,nn),I] = max(PPcormap(n,nn,:),[],3);
            PPmaxlags(n,nn) = lags(I);
        end  
        
        %Negative-Total Waveform
        for nn = 1:Lc
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Nff(:,1,nnn),nf(:,nn,nnn),'coeff'); 
               
               Ncormap(n,nn,:) = squeeze(Ncormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [Nmaxcor(n,nn),I] = max(Ncormap(n,nn,:),[],3);
            Nmaxlags(n,nn) = lags(I);
        end 

        %Negative-Negative Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Nnf = nfsignal;
            Nnf(nfsignal >= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Nff(:,1,nnn),Nnf(:,nnn),'coeff'); 
               
               NNcormap(n,nn,:) = squeeze(NNcormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [NNmaxcor(n,nn),I] = max(NNcormap(n,nn,:),[],3);
            NNmaxlags(n,nn) = lags(I);
        end 
        
        %Positive-Negative Waveform
        for nn = 1:Lc
            nfsignal = squeeze(nf(:,nn,:));
            nfmean = mean(nf(:));
            Pnf = nfsignal;
            Pnf(nfsignal <= nfmean) = nfmean;
            for nnn = 1:NsBS
               %compute nearfield to farfield correlations
               [c,ilags] = xcorr(Nff(:,1,nnn),Pnf(:,nnn),'coeff'); 
               
               PNcormap(n,nn,:) = squeeze(PNcormap(n,nn,:))+c/NsBS;
               lags = ilags/pp.FS;

            end
            [PNmaxcor(n,nn),I] = max(PNcormap(n,nn,:),[],3);
            PNmaxlags(n,nn) = lags(I);
        end 
    end
    
    %Reshape and Sort spatial grid
    [X Y Tmaxcor Tmaxlags Tcormap] = ReshapeGrid(X,Y,Tmaxcor,Tmaxlags,Tcormap);
    [X Y Pmaxcor Pmaxlags Pcormap] = ReshapeGrid(X,Y,Pmaxcor,Pmaxlags,Pcormap);
    [X Y Nmaxcor Nmaxlags Ncormap] = ReshapeGrid(X,Y,Nmaxcor,Nmaxlags,Ncormap);    
    [X Y PPmaxcor PPmaxlags PPcormap] = ReshapeGrid(X,Y,PPmaxcor,PPmaxlags,PPcormap);
    [X Y NNmaxcor NNmaxlags NNcormap] = ReshapeGrid(X,Y,NNmaxcor,NNmaxlags,NNcormap);
    [X Y PNmaxcor PNmaxlags PNcormap] = ReshapeGrid(X,Y,PNmaxcor,PNmaxlags,PNcormap);
end