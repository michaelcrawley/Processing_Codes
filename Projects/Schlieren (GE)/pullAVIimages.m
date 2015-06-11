function pullAVIimages()
    [flist,src_dir] = uigetfile('.avi','Select avi files','Multiselect','on');
    flist = cellfun(@(file) fullfile(src_dir,file),flist,'UniformOutput',0);

    for n = 1:length(flist)
        [~,filename] = fileparts(flist{n});
        mkdir([src_dir,filename]); %make output directory
        
        obj = mmreader(flist{n}); %get mmreader object
%         gimg = zeros(obj.Height,obj.Width,obj.NumberofFrames);
        img = read(obj);
        for nn = 1:obj.NumberOfFrames
            gimg = rgb2gray(img(:,:,:,nn)); %convert to grayscale
            num = num2str([zeros(1,length(num2str(obj.NumberOfFrames))-length(num2str(nn))),nn]);
            imwrite(gimg,[src_dir,filename,filesep,'raw_',num(num~=' '),'.tif'],'TIFF','Compression','lzw');%save img to file directory
        end        
    end
end