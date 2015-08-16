function [out] = TimeResolvedPIV(src,piv_fid,acoustic_fid)

    %DAQ constants
    FS = 4e5;
    BS = 24576;
    NCH = 18;
    lCH = 17; %laser signal from DET210
    NFCH = 6:16; %let's get rid of the first near-field microphone
    alpha = 10;
    DS = 4; %the microphones are low-pass filtered at FS/4, so there is no reason to keep all of that data
    width = 512; %number of points (after downsample) on either side of the laser trigger for correlation

    %Read PIV data first, so we know what blocks to throw out (due to bad
    %images)
    piv = load([src,filesep,piv_fid]);
    bad_vec = piv.badvec_chk;
    fid = fopen([src,filesep,acoustic_fid],'r');
    raw = fread(fid,'float32'); 
    fclose(fid);
    raw = reshape(raw,BS,NCH,[]);
    signal = squeeze(raw(:,lCH,bad_vec));
    NF = raw(:,6:16,piv.badvec_chk); %fix to test how important 1st microphone is
    
    %Identify trigger indices
    trigger = diff(signal);
    trigger1 = trigger/std(trigger(:)) > alpha;
    trigger2 = signal/std(signal(:)) > alpha/4;
    trigger = trigger1.*trigger2(2:end,:);    
    [~,test] = max(trigger);
    test = test + 1;
    
    %Grab relevant signal data
    Nblock = length(test);
    xt = zeros(Nblock,(2*width+1)*size(NF,2));
    for n = 1:Nblock
        tmp = NF(test(n)-DS*width:DS:test(n)+DS*width,:,n);
        

        xt(n,:) = tmp(:)';
    end
    
    %Grab relevant output data
    %FOR NOW, WE ARE JUST GOING TO USE Ux
    d = reshape(piv.data(1).U,[],Nblock)';
    
    %Free up some RAM...
    clear trigger* piv signal NF raw
    
    %Normalize 
    xt_norm = std(xt(:));
    d_norm = max(d(:));
    xt = xt/xt_norm;
    d = d/d_norm;
    
    %ANN
    keyboard
    [nneFun, MSE, weights] = FFN_BP(xt,d,64,'amplitude',1.1);
    
    
    fid = fopen([src,filesep,acoustic_fid],'r');
    raw = fread(fid,'float32'); 
    fclose(fid);
    raw = reshape(raw,BS,NCH,[]);
    NF = raw(:,NFCH,bad_vec);
    NF = NF(:,:,1); %lets just use the first block...
    np = length(NF);
    
    counter = 1;
    for n = DS*width+1:DS:np-DS*width, 
        tmp = NF(n-DS*width:DS:n+DS*width,:)/xt_norm; 
        out(:,counter) = nneFun(tmp(:)')*d_norm; 
        counter = counter+1;
    end
    

    for n = 1:size(out,2)
        pcolor(x,y,out2(:,:,n)); 
        shading interp; colormap jet; colorbar; caxis(clims);
        xlabel('x/D'); ylabel('y/D');
        title(['t = ',num2str(t(n)*1e6),' \mus']);
        
        frame = getframe;
        writeVideo(aviobj,frame);
    end
    close(aviobj);
    
    h = figure;
    movegui(h, 'onscreen');
    rect = get(h,'Position'); 
    rect(1:2) = [0 0]; 
    aviobj = VideoWriter('test.avi');
    open(aviobj);
    for n = 1:size(out2,2)
        pcolor(x,y,out2(:,:,n)); 
        shading interp; colormap jet; colorbar; caxis(clims);
        xlabel('x/D'); ylabel('y/D');
        title(['t = ',num2str(t(n)*1e6),' \mus']);
        
        movegui(h, 'onscreen');
        hold all;

        drawnow;
        writeVideo(aviobj,getframe(gcf,rect));
    end
    close(aviobj);

end



% %%%%%%plot 12/11 channel input results
% for n = 1:size(out1,3)
%     subplot(2,1,1);
%     pcolor(x,y,out1(:,:,n)); shading interp; colormap jet; colorbar; caxis([-50 350]);
%     subplot(2,1,2);
%     pcolor(x,y,out2(:,:,n)); shading interp; colormap jet; colorbar; caxis([-50 350]);
%     drawnow;
% end