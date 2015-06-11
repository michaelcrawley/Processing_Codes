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