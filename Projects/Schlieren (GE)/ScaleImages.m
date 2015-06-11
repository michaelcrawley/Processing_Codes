function ScaleImages(src_dir,flist,backimg,cHi,cWi)

    alpha = 1.1;
    files = cellfun(@(file) fullfile(src_dir,file),flist,'UniformOutput',0);
    
    backimg = double(imread(backimg));
    
    for n = 1:length(files)
        [~,filename] = fileparts(files{n});
        raw = double(imread(files{n})); %read in image, convert data type to double from uint8
        
        avgbck = mean(mean(backimg(500:end,600:end))); %calculate background image intensity scaling parameter
        avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter
        
        sub = raw-alpha*avgimg/avgbck*backimg; %subtract scaled background
        sub(sub < 0) = 0;
        
        minI = min(min(sub(cHi,cWi))); %find minimum illumination
        maxI = max(max(sub(cHi,cWi))); %find maximum illumination
        
        scaled = 255*(sub-minI*ones(size(sub)))/(maxI-minI); %scale image so illumination index goes from 0 to 255
        imwrite(uint8(scaled),[src_dir,filename,'_scaled.tif'],'TIFF','Compression','lzw');        
    end
end