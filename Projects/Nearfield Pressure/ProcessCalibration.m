function ProcessCalibration(src_dir)
    
    %get all Calibration files in directory
    flist = getfiles('CAL*.bin',src_dir);
    
    for n = 1:length(flist)
       [~,fname] = fileparts(flist{n}); %get filename
       md = getMetaFromStr(fname,'MicCalibrationParams'); %get meta data from filename       
       Prms = (2e-5)*10^(md.dB.value/20); %compute rms pressure based off of reference value
       
       %read in data
       fid = fopen([src_dir filesep flist{n}],'r');
       v = fread(fid,'float32');
       fclose(fid);
       
       %Compute rms voltage and conversions
       vrms = sqrt(mean(v.^2));
       PaV = Prms/vrms;
       
       %save output
       save([src_dir filesep fname '.mat'],'vrms','PaV');        
    end

end