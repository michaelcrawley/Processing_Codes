function [Prms X Y] = computeIntensity(src_dir,flist)

    %Temporary constants
    pp.NFCh = 1:8;
    pp.FFCh = 1:11;

    %Initialize variables
    Lf = length(flist);
    Lnfc = length(pp.NFCh);
    Lffc = length(pp.FFCh);
    Prms.nf.filt = zeros(Lf,Lnfc);
    Prms.nf.ufilt = Prms.nf.filt;
    Prms.ff.filt = zeros(Lf,Lffc);
    Prms.ff.ufilt = Prms.ff.filt;
    X = Prms.nf.filt;
    Y = Prms.nf.filt;
    
    h = waitbar(0,'Calculating...0%');
    for n = 1:length(flist)
        waitbar(n/length(flist),h,['Calculating...',num2str(round(100*n/length(flist))),'%']);
        
        %load data
        data = load([src_dir filesep flist{n}]);
        
        %determine x and y locations for microphone signals
        X(n,:) = data.phys.x;
        Y(n,:) = data.phys.y;
        
        %Compute pressure intensities
        Prms.nf.filt(n,:) = mean(sqrt(mean(data.nf.pblocks.smp.^2,1)),3);
        Prms.ff.filt(n,:) = mean(sqrt(mean(data.ff.pblocks.smp.^2,1)),3);   
        Prms.nf.ufilt(n,:) = mean(sqrt(mean(data.nf.pblocks.p.^2,1)),3);
        Prms.ff.ufilt(n,:) = mean(sqrt(mean(data.ff.pblocks.p.^2,1)),3); 
    end
    close(h);
    
    %Reshape and Sort spatial grid
    [X Y Prms.nf.filt Prms.nf.ufilt] = ReshapeGrid(X,Y,Prms.nf.filt,Prms.nf.ufilt);   
end