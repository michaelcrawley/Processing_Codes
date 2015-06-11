function SaveAcousticData(filename,results,data,pp,mat_dir)

    tempstr = strrep(strrep(strrep(strrep(filename(1:end-4),'.','_'),'-','n'),'(','_'),')','_');
    tempvar = genvarname(tempstr);
    tempdata = struct('processing_params',pp,'data',data,'results',results);
    eval([tempvar '=tempdata;']);
    save([mat_dir '\' tempstr '.mat'],tempstr);

end