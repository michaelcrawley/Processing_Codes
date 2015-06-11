function PreProcessNumericalData(flist,src_dir,out_dir)
    %column headers
    headers = {'i','j','k','x','y','z','u','v','w','p','rho'};
    N = length(headers);
    sNf = num2str(length(flist));

    for n = 1:length(flist)
        clc;
        disp(['Processing file: ',num2str(n),' of ',sNf]);
        %open file
        fid = fopen([src_dir filesep flist{n}]);
        [~,fname] = fileparts(flist{n});
        %get iteration time
        theader = fgetl(fid);
        tsplit = regexp(theader,'\s+','split');
        t = str2num(tsplit{2});
        %read data
        raw = cell2mat(textscan(fid,'%n','headerlines',1));
        fclose(fid);
        %reshape data
        for q = 1:N
            eval([headers{q},' = raw(',num2str(q),':N:end);']);
        end
        [~,I,J] = unique(x);
        Lx = length(I);
        Lr = find(J == 1,1,'last');
        for q = 1:N
            eval([headers{q},' = reshape(',headers{q},',Lr,Lx);']);
        end
        %Save data
        save([out_dir filesep fname '.mat'],'t','i','j','k','x','y','z','u','v','w','p','rho');
    end

end