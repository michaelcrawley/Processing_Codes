function [cohmap f X Y] = computeCOH(src_dir,flist)


    %Temporary constants
    pp.BS = 81920;
    pp.sBS = 8192;
    pp.NB = 10;
    pp.NCh = 20;
    pp.NFCh = 13:20; %nearfield microphone channel
    pp.FFCh = 1:11; %farfield microphone channels  
    pp.tCh = 12; %channel number of trigger signal
    pp.CCh = 2; %farfield channel for coherence computation (30 deg)
    pp.FS = 200000; %data acquisition frequency, in Hz
    pp.D = 0.0254;
    pp.AT = 26.1; %atmospheric temperature, in C
    pp.sp = 1;
    pp.nm = 8;
    pp.R = [2.45 2.57 2.69 2.90 3.19 3.46 3.68 3.28 3.15 2.63 2.56]; %radial distance of farfield microphones, in meters   
    
    %initialize variables
    Lf = length(flist);
    Lc = length(pp.NFCh);
    nfrq = ceil((pp.sBS+1)/2);
    X = zeros(Lf,Lc);
    Y = X;
    cohmap = zeros(Lf,Lc,nfrq);
    
    for n = 1:Lf    
        %determine x and y locations for microphone signals
        [~,fname] = fileparts(flist{n});
        phys = NF_ReadFileName(fname,pp);
        X(n,:) = phys.x;
        Y(n,:) = phys.y;
       
        %Read file
        fid = fopen([src_dir filesep flist{n}],'r');
        raw = fread(fid,'float32');fclose(fid);
        data = reshape(raw,pp.BS,pp.NCh,[]); %reshape data into Points x Channels x Blocks
        data = permute(data,[1 3 2]); %reshape data into Points x Blocks x Channels
       
        %compute coherence
        for nn = pp.NFCh
            [cohmap(n,nn-pp.NFCh(1)+1,:) f] = coherence(data(:,:,nn),data(:,:,2),pp.FS,pp.sBS,@hann);          
        end
    end
end