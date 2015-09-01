function PreprocessPIV(src,wfm,reported_dt,ix,iy)
    
    if ~exist('ix','var'), ix = []; end
    if ~exist('iy','var'), iy = []; end

    %Read in LeCroy waveform to determine true dt, for proper velocity
    %scaling
    [laser,info] = readLeCroyWaveJetWfm(wfm);
    [~,I] = max(sum(laser));
    laser = laser(:,I); %I am assuming that nothing else is plugged into the WaveJet aside from the DET210
    t = (0:length(laser)-1)*info.dx;
    [vals,locs] = findpeaks(laser,'minpeakdistance',2,'sortstr','descend','npeaks',2);
    figure;
    plot(t*1e6,laser,t(locs)*1e6,vals,'*'); %convert time to us
    xlabel('time (\mus)');
    true_dt = -diff(t(locs))*1e6;
    disp(['Calculated true dt = ',num2str(true_dt),' (micro-seconds)']);
    chk_flag = input('Is this correct? [y/n]: ','s');
    if ~(strcmpi(chk_flag,'y') || strcmpi(chk_flag,'yes'))
        true_dt = input('Enter correct dt in micro-seconds: ');
    end
    scaling = reported_dt/true_dt;
    
    %Read in all images, throw out bad ones
    [flist,~,folders] = getfiles('data.mat',src,'-a');
    
    keyboard;
   
    
    for n = 1:length(folders)
        
        %Load data
        load([src,filesep,flist{n}{1}]);

        for q = 1:length(data)
            %Identify bad images based on number of vectors
            nvec = squeeze(sum(sum((data(q).U+data(q).V > 0),2),1));
            chk1 = (nvec/numel(data(q).X)) > 0.75;
            
            %Identify bad images based on jet core velocity averages
            test = zeros(size(data(q).U,3),1);
            for k = 1:size(data(1).U,3)
                tmp = data(1).U(ix,iy,k);
                [~,test(k,1)] = nzstats(tmp(:),1);
            end
            chk2 = test < (mean(test)+1.25*std(test)) & test > (mean(test)-1.25*std(test));
            badvec_chk(:,q) = chk1 & chk2;            
            
            data(q).U = data(q).U(:,:,badvec_chk(:,q));
            data(q).V = data(q).V(:,:,badvec_chk(:,q));  
            
            %Rescale per correct laser dt
            data(q).U = scaling*data(q).U;
            data(q).V = scaling*data(q).V; 
        end

        %Save in data directory
        save([folders{n},filesep,'preprocessed.mat'],'data','badvec_chk','true_dt','scaling');
    end
    
end