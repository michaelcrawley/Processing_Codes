function TimeResolvedPIV_Velocity2(src,piv_fid,acoustic_fid,arch,tag)
    outdir = 'Test Reconstructions';

    %DAQ constants
    FS = 4e5;
    BS = 24576;
    NCH = 18;
    lCH = 17; %laser signal from DET210
    NFCH = 5:16;
    alpha = 10;
    DS = 2; %the microphones are low-pass filtered at FS/4, so there is no reason to keep all of that data
    width = 512; %number of points (after downsample) on either side of the laser trigger for correlation
    Nblocks = 1500;
    Amplitude = 1.1;

    %Read PIV data first, so we know what blocks to throw out (due to bad
    %images)
    piv = load([src,filesep,piv_fid]);
    fid = fopen([src,filesep,acoustic_fid],'r');
    raw = fread(fid,'float32'); 
    fclose(fid);
    raw = reshape(raw,BS,NCH,[]);  
    Nimag = size(raw,3);
    missed_blocks = Nblocks-Nimag;
    signal = squeeze(raw(:,lCH,piv.badvec_chk(missed_blocks+1:end,1)));
    NF = raw(:,NFCH,piv.badvec_chk(missed_blocks+1:end,1));
    
    %find images corresponding to missed blocks
    cut = sum(ones(missed_blocks,1) .* piv.badvec_chk(1:missed_blocks,1));
    piv.data(1).U = piv.data(1).U(:,:,1+cut:end);
    piv.data(1).V = piv.data(1).V(:,:,1+cut:end);
    

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
    phi_chk = piv.data(1).Y(1,:,1) > 0;
    U = piv.data(1).U(:,phi_chk,:);
    V = piv.data(1).V(:,phi_chk,:);
    
    [~,Um] = nzstats(U,3);
    [~,Vm] = nzstats(V,3);
    
    Ufluct = U - repmat(Um,[1 1 Nblock]);
    Vfluct = V - repmat(Vm,[1 1 Nblock]);
%     uf_norm = max(abs(Ufluct(:)));
%     vf_norm = max(abs(Vfluct(:)));
    
    d = [reshape(Ufluct,[],Nblock)',reshape(Vfluct,[],Nblock)'];
    
    %Free up some RAM...
    x = piv.data(1).X(:,phi_chk);
    y = piv.data(1).Y(:,phi_chk);
    inputs = NF(:,:,1);
    clear trigger* piv signal raw U V NF Ufluct Vfluct
    
    %Normalize 
    xt_norm = max(abs(xt),[],2);
    d_norm = max(abs(d),[],2);
    xt = xt/xt_norm;
    d = d/d_norm;
    
    %ANN
    [nne, mse, weights] = FFN_BP(xt,d,arch,'amplitude',Amplitude,'maxepoch',5e3,'lrp',0.002);
    
    fname = ['FFNBP_arch',num2str(arch),tag,'UVf.mat'];    
    save([outdir filesep fname],'-v7.3','nne','mse','BS','DS','width','d_norm','xt_norm','d','xt','inputs','arch','Um','Vm','x','y','src','acoustic_fid','piv_fid','weights','Amplitude');

end