function LES_PreprocessPhavg(src,out)

    %Constant fixes
    ngheaders = 6; %there are six columns in the grid file: i, j, k, x, y, z
    ndheaders = 5; %there are five columns in the data file: u, v, w, p, rho
    
    %%%%    Grid File    
    %%%%%%%%%%%%%%%%%
    
    %Read in grid file data
    grd = getfiles('*.grd',src);
    fid = fopen([src filesep grd{1}],'r');
    headers = textscan(fid,'%s',ngheaders); 
    data = textscan(fid,repmat('%f ',1,ngheaders));    
    fclose(fid);
    
    %parse & reshape
    np = cellfun(@(x) length(unique(x)),data(1:ngheaders/2));%first half of headers are indices, second half are physical locations    
    [~,order] = sort(cellfun(@(x) find(diff(x),1),data(1:ngheaders/2))); 
    np = np(order);
    disp('Processing grid file...');
    for n = 1:ngheaders  
        grid.(headers{1}{n}) = permute(reshape(data{n},np),order);
        if n <= ngheaders/2
            grid.(headers{1}{n}) = uint16(grid.(headers{1}{n})); %lets save some space
        end
    end
    [grid.theta,grid.r] = cart2pol(grid.y,grid.z); %get polar coordinates as well   
    grid.theta(1,1,:) = grid.theta(2,2,:); %fix so centerline isn't pure zero
    
    save([out,filesep,'grid.mat'],'grid');
    clear data 

    %%%%    Data Files
    %%%%%%%%%%%%%%%%%%
    disp('Processing Phase-averaged files...');
    
    %Find files, sort on time index    
    flist = getfiles('*.txt',src);
    for n = 1:length(flist)
        [~,fname] = fileparts(flist{n});
        
        fid = fopen([src filesep flist{n}],'r');
        %get time
        tmp = textscan(fid,'%s %f',1);
        t = tmp{2};
        
        headers = textscan(fid,'%s',ndheaders); 
        tmp = textscan(fid,repmat('%f ',1,ngheaders));    
        fclose(fid);
        
        for q = 1:ndheaders  
            data.(headers{1}{q}) = permute(reshape(tmp{q},np),order);
        end
        
        %compute cylindrical components
%         data.u_theta = cos(grid.theta).*data.w + sin(grid.theta).*data.v;
%         data.u_r = -sin(grid.theta).*data.w + cos(grid.theta).*data.v;
        
        data.u_r = sin(grid.theta).*data.w + cos(grid.theta).*data.v;
        data.u_theta = -data.w.*cos(grid.theta) + data.v.*sin(grid.theta);
        
        save([out,filesep,fname,'.mat'],'-struct','data');
    end
end