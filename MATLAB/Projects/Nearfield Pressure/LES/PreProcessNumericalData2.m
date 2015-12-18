function PreProcessNumericalData2(src_dir,out_dir)
%This code will grab the pressure out of the LES files for all time, and
%save them in a way to match the output of NF_Preprocessing.m (i.e. nf
%pressure along a single radial distance for all x and t.

    %%Grab grid info 
    fgrid = getfiles('*.grd',src_dir);
    fid = fopen([src_dir filesep fgrid{1}],'r');
    gridheader = {'i','j','k','x','y','z'};
    N = length(gridheader);
    %Grab locations
    raw = cell2mat(textscan(fid,'%n','headerlines',1)); %header is already read, so we don't need to do it again
    fclose(fid);
    %Reshape grid
    for q = 1:N
        eval([gridheader{q},' = raw(',num2str(q),':N:end);']);
    end
    [~,I,J] = unique(x);
    Lx = length(I);
    Lr = find(J == 1,1,'last');
    for q = 1:N
        eval([gridheader{q},' = reshape(',gridheader{q},',Lr,Lx);']);
    end
    %Transform from y,z to r
    r = sqrt(y.^2+z.^2);
    %Transform to inches
    D = 0.0254;
    x = x/D;
    r = r/D;
    
    %%Read in data
    headers = {'u','v','w','p','rho'};
    N = length(headers);
    fdata = getfiles('*.txt',src_dir);
    Lt = length(fdata);
    t = zeros(Lt,1);
    pressure = zeros(Lr,Lx,Lt);
    h = waitbar(0,'Grabbing Data...');
    for n = 1:Lt
        fid = fopen([src_dir filesep fdata{n}],'r');
        
        %get iteration time
        theader = fgetl(fid);
        tsplit = regexp(theader,'\s+','split');
        t(n) = str2num(tsplit{2});
        
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
        
        %grab pressure data
        pressure(:,:,n) = p;
        
        %update user
        waitbar(n/Lt,h);
    end
    
    %Sort for time index
    [t,I] = sort(t);
    pressure = pressure(:,:,I);
    
    %Save data to individual mat files
    phys.dt = mean(diff(t));
    phys.t = t;
    phys.x = x(1,:);
    waitbar(0,h,'Saving Data...');
    for n = 1:Lr
        nf.LES.p = permute(squeeze(pressure(n,:,:)),[2 1]);
        phys.y = r(n,:);
        save([out_dir,filesep,'j',num2str(n),'.mat'],'nf','phys');
        waitbar(n/Lr,h);
    end
    close(h);


end