function[R] = CAvgImgs(src_dir,flist,backimg,cHi,cWi,dthreshold)

    files = cellfun(@(file) fullfile(src_dir,file),flist,'UniformOutput',0);
    
    backimg = double(imread(backimg));
    avgbck = mean(mean(backimg(500:end,600:end))); %calculate background image intensity scaling parameter
    alpha = 1.1;
    
    M = length(cHi);
    N = length(cWi);
    nF = length(flist);
      
    %%create mean_sub_cropped image
    mean_sub_cropped = zeros(M,N);
    for n = 1:nF
       raw = double(imread(files{n})); 
       avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter
       
       sub = raw-alpha*avgimg/avgbck*backimg; %subtract scaled background
       
       sub_cropped = sub(cHi,cWi); %crop image down to section we are only interested in
       
       mean_sub_cropped = mean_sub_cropped+(1/nF)*sub_cropped;
    end
    
    %Creates matrix for correlations
    chk = zeros(M*N,length(flist)); %initial matrix for correlation coeficients
    for n = 1:nF
       raw = double(imread(files{n})); 
       avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter       
       sub = raw-alpha*avgimg/avgbck*backimg; %subtract scaled background       
       cropped = sub(cHi,cWi); %crop image down to section we are only interested in       
       shifted_cropped = cropped-mean_sub_cropped;       
       chk(:,n) = reshape(shifted_cropped,M*N,1);
    end

    R = corrcoef(chk); %correlation coefficient
     
    %%Find correlated images for phase 0
    done = false; 
    LL = 0.2;
    threshold = dthreshold;
    nImg = nF/20; %define required number of images for correlation as at least 10% of set
    while ~done
        alpha = R > threshold;
        L = sum(alpha);   %Finds the number of images with correlation above "threshold"
        if max(L) > nImg   %If the necessary number of images is found
            done = true;
            CorrCoeff = threshold;
            disp([num2str(max(L)) ' Positive Match Found at Corr. Coeff. = ' num2str(threshold)])
        elseif threshold < LL   %If threshold has fallen below acceptable lower limit
            done = true;
            CorrCoeff = 0;
            disp('No Positive Result')
        else
            threshold = threshold-0.01;
        end
    end
    
    [~,I] = max(L);
    matches0 = alpha(:,I);
    Ph0 = flist(matches0);
    
    %%Find correlated images for phase pi
    done = false;
    threshold = -dthreshold;
    while ~done
        alpha = R(:,I) < threshold;
        L = sum(alpha);
        if L > nImg
            done = true;
            CorrCoeff(2) = threshold;
            disp([num2str(L) ' Negative Match Found at Corr. Coeff. = ' num2str(threshold)])
        elseif threshold > -LL
            done = true;
            CorrCoeff(2) = 0;
            disp('No Positive Result')
        else
            threshold = threshold + 0.01;
        end
    end
    
    matches180 = alpha;
    Ph180 = flist(matches180);
    
    %Creates Scaled, averaged image for phase 0
    CAvgPh0 = zeros(size(backimg));
    for n = 1:length(Ph0)
        raw = double(imread(Ph0{n})); %read in image, convert data type to double from uint8
        avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter
        
        sub = raw-1.1*avgimg/avgbck*backimg; %subtract scaled background
        sub(sub < 0) = 0;
        
        minI = min(min(sub(cHi,cWi))); %find minimum illumination
        maxI = max(max(sub(cHi,cWi))); %find maximum illumination
        
        scaled = 255*(sub-minI*ones(size(sub)))/(maxI-minI); %scale image so illumination index goes from 0 to 255
        CAvgPh0 = CAvgPh0+(1/length(Ph0))*scaled;
    end
    imwrite(uint8(CAvgPh0),[src_dir,filesep,'CAvgPh0.tif'],'TIFF','Compression','lzw'); 
    
    %Creates Scaled, averaged image for phase 180
    CAvgPh180 = zeros(size(backimg));
    for n = 1:length(Ph180)
        raw = double(imread(Ph180{n})); %read in image, convert data type to double from uint8
        avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter
        
        sub = raw-1.1*avgimg/avgbck*backimg; %subtract scaled background
        sub(sub < 0) = 0;
        
        minI = min(min(sub(cHi,cWi))); %find minimum illumination
        maxI = max(max(sub(cHi,cWi))); %find maximum illumination
        
        scaled = 255*(sub-minI*ones(size(sub)))/(maxI-minI); %scale image so illumination index goes from 0 to 255
        CAvgPh180 = CAvgPh180+(1/length(Ph180))*scaled;
    end
    imwrite(uint8(CAvgPh180),[src_dir,filesep,'CAvgPh180.tif'],'TIFF','Compression','lzw');
    
    
    %Calculates autocorrelation for phase 0
    CAvgPh0_auto = zeros(1,2*N-1);
    for n = 1:length(Ph0)
        raw = double(imread(Ph0{n})); %read in image, convert data type to double from uint8
        avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter       
        sub = raw-1.1*avgimg/avgbck*backimg; %subtract scaled background       
        cropped = sub(cHi,cWi); %crop image down to section we are only interested in       
        shifted_cropped = cropped-mean_sub_cropped;  
        CAvgPh0_auto = CAvgPh0_auto + (1/length(Ph0))*xcorr2_1d(shifted_cropped,shifted_cropped,2);
    end
    figure; plot([-N+1:0 1:N-1],CAvgPh0_auto); title('CAvgPh0 autocorrelation'); saveas(gcf,'CAvgPh0 autocorrelation','fig');
        
    %Calculates autocorrelation for phase 180
    CAvgPh180_auto = zeros(1,2*N-1);
    for n = 1:length(Ph180)
        raw = double(imread(Ph180{n})); %read in image, convert data type to double from uint8
        avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter       
        sub = raw-1.1*avgimg/avgbck*backimg; %subtract scaled background       
        cropped = sub(cHi,cWi); %crop image down to section we are only interested in       
        shifted_cropped = cropped-mean_sub_cropped;  
        CAvgPh180_auto = CAvgPh0_auto + (1/length(Ph0))*xcorr2_1d(shifted_cropped,shifted_cropped,2);
    end
    figure; plot([-N+1:0 1:N-1],CAvgPh180_auto); title('CAvgPh180 autocorrelation'); saveas(gcf,'CAvgPh180 autocorrelation','fig');
    close all;

end