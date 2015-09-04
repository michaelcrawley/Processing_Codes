function TimeResolvedPIV_Velocity(src,piv_fid,acoustic_fid,arch,tag)
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

    %Read PIV data first, so we know what blocks to throw out (due to bad
    %images)
    piv = load([src,filesep,piv_fid]);
    fid = fopen([src,filesep,acoustic_fid],'r');
    raw = fread(fid,'float32'); 
    fclose(fid);
    raw = reshape(raw,BS,NCH,[]);    
    signal = squeeze(raw(:,lCH,piv.badvec_chk(:,1)));
    NF = raw(:,NFCH,piv.badvec_chk(:,1));
%     NF = raw(:,6:16,piv.badvec_chk); %fix to test how important 1st microphone is
    
    

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
    xt_norm = max(abs(xt(:)));
    d_norm = max(abs(d(:)));
    xt = xt/xt_norm;
    d = d/d_norm;
    
    %ANN
    [nne, mse, weights] = FFN_BP(xt,d,arch,'amplitude',1.1,'maxepoch',5e3,'lrp',0.002);
    
    fname = ['FFNBP_arch',num2str(arch),tag,'UVf.mat'];    
    save([outdir filesep fname],'-v7.3','nne','mse','BS','DS','width','d_norm','xt_norm','d','xt','inputs','arch','Um','Vm','x','y','src','acoustic_fid','piv_fid','weights');
    
%     counter = 1;
%     for n = DS*width+1:DS:BS-width*DS
%         tmp = inputs(n-DS*width:DS:n+DS*width,:);
%         tmp = tmp(:)'/xt_norm;
%         out(:,counter) = nne(tmp)*d_norm;
%         counter = counter + 1;
%     end
%     
% 
%     for n = 1:size(out2,2)
%         pcolor(x,y,out2(:,:,n)); 
%         shading interp; colormap jet; colorbar; caxis(clims);
%         xlabel('x/D'); ylabel('y/D');
%         title(['t = ',num2str(t(n)*1e6),' \mus']);
%         
%         frame = getframe;
%         writeVideo(aviobj,frame);
%     end
%     close(aviobj);
%     
%     h = figure;
%     movegui(h, 'onscreen');
%     rect = get(h,'Position'); 
%     rect(1:2) = [0 0]; 
%     aviobj = VideoWriter('test.avi');
%     open(aviobj);
%     for n = 1:size(out2,2)
%         pcolor(x,y,out2(:,:,n)); 
%         shading interp; colormap jet; colorbar; caxis(clims);
%         xlabel('x/D'); ylabel('y/D');
%         title(['t = ',num2str(t(n)*1e6),' \mus']);
%         
%         movegui(h, 'onscreen');
%         hold all;
% 
%         drawnow;
%         writeVideo(aviobj,getframe(gcf,rect));
%     end
%     close(aviobj);

end