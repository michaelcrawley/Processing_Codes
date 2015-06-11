function [tstampnum,tstampstr] = gettimestamp(src_dir,flist)
%gets timestamp data for specified files.  tstampnum is the serial number,
%whereas tstampstr is the formatted string output.

    tmp = dir(src_dir);
    allfiles = {tmp.name};
    dates = {tmp.date};
    filematch = cellfun(@(file) strmatch(file,flist),allfiles,'UniformOutput',0);
    
    for n = 1:length(filematch)
       if isempty(filematch{n}), filematch{n} = 0; end      
    end
    filematch = cell2mat(filematch);
    
    tstampstr = cell(1,length(flist));
    for n = 1:length(flist)
        [~,I] = find(filematch == n,1);
        tstampstr{n} = dates{I};
    end
    tstampnum = datenum(tstampstr);
end