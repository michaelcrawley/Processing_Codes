function TimeResolvedPIV_Pressure(src,piv_fid,acoustic_fid,arch,tag)
    outdir = 'Test Reconstructions';

    %DAQ constants
    FS = 4e5;
    BS = 24576;
    NCH = 18;
    lCH = 17; %laser signal from DET210
    NFCH = 5:16;
    alpha = 10;
    DS = 4; %the microphones are low-pass filtered at FS/4, so there is no reason to keep all of that data
    width = 256; %number of points (after downsample) on either side of the laser trigger for correlation

    %Read PIV data first, so we know what blocks to throw out (due to bad
    %images)
    piv = load([src,filesep,piv_fid]);
    fid = fopen([src,filesep,acoustic_fid],'r');
    raw = fread(fid,'float32'); 
    fclose(fid);
    raw = reshape(raw,BS,NCH,[]);    
    signal = squeeze(raw(:,lCH,piv.badvec_chk(:,1)));
    NF = raw(:,NFCH,piv.badvec_chk(:,1));

    %Identify trigger indices
    trigger = diff(signal);
    trigger1 = trigger/std(trigger(:)) > alpha;
    trigger2 = signal/std(signal(:)) > alpha/4;
    trigger = trigger1.*trigger2(2:end,:);    
    [~,test] = max(trigger);
    test = test + 1;
    
    %Grab relevant signal data
    Nblock = length(test);
    xt = zeros(Nblock,(2*width+1)*length(NFCH));
    for n = 1:Nblock
        tmp = NF(test(n)-DS*width:DS:test(n)+DS*width,:,n);
        xt(n,:) = tmp(:)';
    end
    
    %Grab relevant output data
    [M,N,L] = size(piv.sol_P);
ex    spatial_cutout = repmat(tukeywin(N,.05).',[M 1]).*exp(-(2e5)*(piv.r.^4));
    for n = 1:L
        piv.sol_P(:,:,n) = piv.sol_P(:,:,n).*spatial_cutout;
    end        
    d = reshape(piv.sol_P,[],Nblock)';
    
    %Free up some RAM...
    r = piv.r;
    z = piv.z;
    inputs = NF(:,:,1);
    clear trigger* piv signal raw U V NF Ufluct Vfluct
    
    %Normalize 
    xt_norm = std(xt(:));
    d_norm = std(d(:));
    xt = xt/xt_norm;
    d = d/d_norm;
    
    %ANN
    [nne, mse, weights] = FFN_BP(xt,d,arch,'amplitude',1.1,'maxepoch',5e3,'lrp',0.002);
    
    fname = ['FFNBP_arch',num2str(arch),tag,'_Sol_P.mat'];    
    save([outdir filesep fname],'-v7.3','nne','mse','BS','DS','width','d_norm','xt_norm','d','xt','inputs','z','r','src','acoustic_fid','piv_fid','weights');

end