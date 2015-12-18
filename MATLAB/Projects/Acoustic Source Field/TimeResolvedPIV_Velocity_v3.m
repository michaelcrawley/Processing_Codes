function TimeResolvedPIV_Velocity_v3(src,piv_fid,acoustic_fid,arch,tag,mask)
    outdir = 'Test Reconstructions';
    mkdir(outdir);
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
    matversion = 4;
    nmodes = 1000;
    if ~exist('mask','var'), mask = 1; end

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
    U = piv.data(1).U;
    V = piv.data(1).V;    
    [~,Um] = nzstats(U,3);
    [~,Vm] = nzstats(V,3);
    x = piv.data(1).X;
    y = piv.data(1).Y;    
    R = zeros(numel(x)*2,Nblock); 
    inputs = NF(:,:,1);
    Um = Um.*mask';
    Vm = Vm.*mask';
    
    for n = 1:Nblock
        ut = U(:,:,n).*mask' - Um;
        vt = V(:,:,n).*mask' - Vm;
        R(:,n) = [ut(:);vt(:)];
    end       
    clear trigger* piv signal raw NF U V Ufluct Vfluct%Free up some RAM...
    
    [phi,lambda,ak] = SnapShotPOD(R,false);
    phi = phi(:,1:nmodes);
    ak = ak(1:nmodes,:);

    %Normalize 
    xt_norm = max(abs(xt(:)));
    xt = xt/xt_norm;  
    d_norm = max(abs(ak(:)));
    d = ak'/d_norm;
        
    %Save unnecessary outputs before iterating
    fname = ['FFNBP_arch',num2str(arch),tag,'UVf_POD.mat'];    
    save([outdir filesep fname],'-v7.3','BS','DS','width','d_norm','xt_norm','d','xt','inputs','arch','Um','Vm','x','y','src','acoustic_fid','piv_fid','matversion','phi','ak','lambda','nmodes');
    clear BS DS width d_norm xt_norm inputs Um Vm x y phi ak lambda R;
    
    %ANN
    [nne, mse, weights] = FFN_BP(xt,d,arch,'amplitude',1.05,'maxepoch',1e4,'lrp',0.02);    
    save([outdir filesep fname],'-append','nne','mse','weights');

end