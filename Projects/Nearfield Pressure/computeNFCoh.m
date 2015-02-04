function [cohmap f X Y phase cpsd xpsd] = computeNFCoh(src_dir,flist,save_name)
%This function is designed to work on the filtered (non-averaged)
%microphone data. It is designed to work with _proc mat-files 1.1.

    %Processing channel constants
    pp.NFCh = 1:8; %nearfield microphone channel for correlations
    
    %DAQ constants
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.sBS = 8192; %processing block size
    
    %initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh)-1;
    nfrq = ceil((pp.sBS+1)/2);
    X = zeros(Lf,Lc);
    Y = X;
    cohmap = zeros(Lf,Lc,nfrq);
    phase = cohmap;
    cpsd = cohmap;
    xpsd = cohmap;
    
    h = waitbar(0,'Calculating...0% Complete');
    for n = 1:Lf  %file number  
        %determine x and y locations for microphone signals
        data = load([src_dir filesep flist{n}]);
        X(n,:) = data.phys.x(1:end-1);
        Y(n,:) = data.phys.y(1:end-1);
        
        for nn = 1:Lc %coh channel number
            x1 = squeeze(data.nf.pblocks.smp(:,nn,:));
            x2 = squeeze(data.nf.pblocks.smp(:,nn+1,:));
            [cohmap(n,nn,:), f, phase(n,nn,:), cpsd(n,nn,:), xpsd(n,nn,:)] = coherence(x2,x1,pp.FS,pp.sBS,@hann);               
        end
        waitbar(n/Lf,h,['Calculating...',num2str(round(n/Lf*100)),'% Complete']);
    end
    close(h);
    
    %Reshape matrices
    [X Y cohmap phase cpsd xpsd] = ReshapeGrid(X,Y,cohmap,phase,cpsd, xpsd);
    
    %Save Data
    save(save_name,'X','Y','cohmap','phase','cpsd','xpsd','f','pp');
end