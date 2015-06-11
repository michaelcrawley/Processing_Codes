function[R] = CondAverageImg(src_dir,flist,cHi,cWi,backimg,Rthreshold)
%Inputs:
%   src_dir:    source directory for images 
%   flist:      list of files for conditional-averaging (cell array)
%   cHi & cWi:  indices for correlation window (Height and Width)
%   backimg:    background image (full file name)
%   Rthreshold: starting correlation coefficient threshold (optional)

    if ~exist('Rthreshold','var'), Rthreshold = 0.9; end
       
    %Scale Raw Images
    scaled_dir = [src_dir filesep 'scaled' filesep];
    mkdir(scaled_dir);
    ScaleImages(src_dir,flist,backimg,cHi,cWi,scaled_dir,1.1);
    
    %%Find Image Matches based on correlation coefficient
    Cavg_dir = [src_dir filesep 'Cavg' filesep];
    mkdir(Cavg_dir);
    scaled_flist = PullDirFiles('*.tif',scaled_dir);
    
    %Find Image Pairs
    M = length(cHi);
    N = length(cWi);
    nF = length(flist);
    
    chk = zeros(M*N,length(scaled_flist));
    for n = 1:nF
        img = double(imread([scaled_dir filesep scaled_flist{n}]));
        chk(:,n) = reshape(img(cHi,cWi),M*N,1);
    end
    
    %%Binning
    chk = round(chk/2);
    chk(chk > 128) = 128;
%     chk = double(chk >= 0.5*255);
    
    R = corrcoef(chk); %correlation coefficient
    
    done = false; 
    LL = 0.2;
    threshold = Rthreshold;
    nImg = nF/10; %define required number of images for correlation as at least 10% of set
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
    Ph0 = scaled_flist(matches0);
    
    %Save Correlated images to Cavg directory
    cellfun(@(file) copyfile([scaled_dir filesep file],[Cavg_dir filesep file]),Ph0);
    AverageImages(Cavg_dir,Ph0,['CAvg Ph0,alpha',num2str(round(threshold*100))]);

end

function ScaleImages(src_dir,flist,backimg,cHi,cWi,out_dir,alpha)

    if ~exist('alpha','var'), alpha = 1.0; end 
    
    backimg = double(imread(backimg));
    
    for n = 1:length(flist)
        [~,filename] = fileparts(flist{n});
        filenum = filename(regexp(filename,'_')+1:end);
        raw = double(imread([src_dir filesep flist{n}])); %read in image, convert data type to double from uint8
        
        avgbck = mean(mean(backimg(500:end,600:end))); %calculate background image intensity scaling parameter
        avgimg = mean(mean(raw(500:end,600:end))); %calculate raw image intensity scaling parameter
        
        sub = raw-alpha*avgimg/avgbck*backimg; %subtract scaled background
        sub(sub < 0) = 0;
        
        minI = min(min(sub(cHi,cWi))); %find minimum illumination
        maxI = max(max(sub(cHi,cWi))); %find maximum illumination
        
        scaled = 255*(sub-minI*ones(size(sub)))/(maxI-minI); %scale image so illumination index goes from 0 to 255
        imwrite(uint8(scaled),[out_dir,'scaled_',filenum,'.tif'],'TIFF','Compression','lzw');        
    end
end

function [flist,src_dir] = PullDirFiles(condition,src_dir)
    if ~exist('condition','var'), condition = ''; end    
    if ~exist('src_dir','var'), src_dir = uigetdir('Specify Directory'); end
    
    list = dir([src_dir,filesep,condition]);
    chk = [list.isdir];
    flist = {list(~chk).name};    
end

function AverageImages(src_dir,flist,savename)

    flist = cellfun(@(file) fullfile(src_dir,file),flist,'UniformOutput',0);
    info = imfinfo(flist{1});
    [~,savename] = fileparts(savename);
      
    avgimg = zeros(info.Height,info.Width);
    for n = 1:length(flist)
       avgimg = avgimg+(1/length(flist))*double(imread(flist{n},'TIFF')); 
    end
    
    avgimg = uint8(avgimg);
    imwrite(avgimg,[src_dir,filesep,savename,'.tif'],'TIFF','Compression','lzw');    
end