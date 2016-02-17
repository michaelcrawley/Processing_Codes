function voltage = Read_LoadCell(fname,version)

    fid = fopen(fname);
    if version == 8
        tmp = textscan(fid,'%f','headerlines',1);
        voltage = reshape(tmp{1},6,[]);
    else
        tmp = fread(fid,'float32');
        voltage = reshape(tmp,[],6)';
    end
    fclose(fid);
    
end