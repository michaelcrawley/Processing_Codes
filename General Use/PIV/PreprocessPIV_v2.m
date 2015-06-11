function PreprocessPIV(src,wfm,reported_dt)
    

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
    for n = 1:length(folders)
        
        %Load data
        tmp = load([src filesep flist{n}{1}]);

        %Identify bad images based on number of vectors
        nvec = squeeze(sum(sum((tmp.data(1).U+tmp.data(1).V > 0),2),1));
        badvec_chk = (nvec/numel(tmp.data(1).X)) > 0.75;
        data.U = tmp.data(1).U(:,:,badvec_chk);
        data.V = tmp.data(1).V(:,:,badvec_chk);  

        %Rescale per correct laser dt
        data.U = scaling*data.U;
        data.V = scaling*data.V; 
        
        data.X = tmp.data(1).X;
        data.Y = tmp.data(1).Y;

        %Save in data directory
        save([folders{n},filesep,'preprocessed.mat'],'data','badvec_chk','true_dt','scaling');
    end
    
=======
function PreprocessPIV(src,wfm,reported_dt)
    

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
    for n = 1:length(folders)
        
        %Load data
        load([flist{n}{1}]);

        for q = 1:length(data)
            %Identify bad images based on number of vectors
            nvec = squeeze(sum(sum((data(q).U+data(q).V > 0),2),1));
            badvec_chk = (nvec/numel(data(q).X)) > 0.75;
            data(q).U = data(q).U(:,:,badvec_chk);
            data(q).V = data(q).V(:,:,badvec_chk);  
            
            %Rescale per correct laser dt
            data(q).U = scaling*data(q).U;
            data(q).V = scaling*data(q).V; 
        end

        %Save in data directory
        save([folders{n},filesep,'preprocessed.mat'],'data','badvec_chk','true_dt','scaling');
    end
    
>>>>>>> cd808df814a9f73427ee7fbf4f601c1900a034ce
end