function LES_PreprocessSnapshots(src,out,dsfactor)
%Note that the standard definition, per Rachelle's LES simulation mesh, is
%that the jet axis dimension is labeled as 'x', not 'z'
    %Constant fixes
    folders = regexp(src,filesep,'split');
    local = folders{end};
    switch lower(local(1))
        case 'k'
            order = [2 1 3];
            ngheaders = 6; %there are six columns in the grid file: i, j, k, x, y, z
            ndheaders = 5; %there are five columns in the data file: u, v, w, p, rho
    end    
    
    %%%%    Grid File    
    %%%%%%%%%%%%%%%%%
    
    %Read in grid file data
    grd = getfiles('*.grd',src,'-a');
    fid = fopen([src filesep grd{1}],'r');
    headers = textscan(fid,'%s',ngheaders); 
    data = textscan(fid,repmat('%f ',1,ngheaders));    
    fclose(fid);
    
%     %parse & reshape
%     np = cellfun(@(x) length(unique(x)),data(1:ngheaders/2));%first half of headers are indices, second half are physical locations    
%     np = np(order);
%     disp('Processing grid file...');
%     for n = 1:ngheaders  
%         grid.(headers{1}{n}) = reshape(data{n},np);
%     end
%     [grid.theta,grid.r] = cart2pol(grid.y,grid.z); %get polar coordinates as well
%     clear data;
    
    %parse & reshape
    np = cellfun(@(x) length(unique(x)),data(1:ngheaders/2));%first half of headers are indices, second half are physical locations   
    tmp = cellfun(@(x) find(diff(x),1),data(1:ngheaders/2),'uniformoutput',false);
    chk = cellfun(@isempty,tmp);
    tmp{chk} = 0;
    [~,order] = sort(cell2mat(tmp));
%     [~,order] = sort(cellfun(@(x) find(diff(x),1),data(1:ngheaders/2))); 
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
    clear data
    
    %%%%    Actuator File
    %%%%%%%%%%%%%%%%%%%%%
    
    %Find actuator file, if it exists
    act = getfiles('*.act',src);
    if ~isempty(act)
        disp('Actuator signal identified...');
        flags.act = true;
        
        %read in data
        fid = fopen([src filesep act{1}],'r');
        signal = textscan(fid,'%f %f');
        fclose(fid);
        trigger = signal{2};      
    else
        warning('No actuator signal found!');
        flags.act = false;
        trigger = [];
    end
    
    %%%%    Data Files
    %%%%%%%%%%%%%%%%%%
    disp('Sorting and downsampling files...');
    
    %Find files, sort on time index    
    flist = getfiles('*.txt',src);
    tmp = regexp(flist,'mv(\d*).txt','tokens','once');
    it = cellfun(@(x) str2double(x),tmp);
    [it,I] = sort(it);
    flist = flist(I);
    if length(unique(diff(it))) > 1 %make sure the files represent continuous samples
        warning('Missing Files Detected!');
        flags.continuous = false;
        chk = input('Do you want to continue? [y/n] ','s');
        if strcmpi(chk,'n') || strcmpi(chk,'no'), return; end  
    else
        flags.continuous = true;
    end
    
    %downsample data in time
    flist = flist(1:dsfactor:end);
    Nf = length(flist);
    if ~isempty(trigger), trigger = trigger(1:dsfactor:end); end
    
    %get header info
    fid = fopen([src filesep flist{1}],'r');
    fgetl(fid); %skip time info        
    headers = textscan(fid,'%s',ndheaders);%get data header line
    fclose(fid);
    
    %%%% Process variables
    %%%%%%%%%%%%%%%%%%%%%%
    disp('Processing variables...');
    cpb = ConsoleProgressBar();
    t = zeros(Nf,1);   
    
%     u_theta = zeros([np(order),Nf]);
%     u_r = zeros([np(order),Nf]);
    
%     header_v = strmatch('v',headers{1});
%     header_w = strmatch('w',headers{1});
    
    for q = 1:ndheaders
        data = zeros([np(order),Nf]);
        cpb.start();
        for n = 1:Nf
            fid = fopen([src filesep flist{n}],'r');

            %get time
            tmp = textscan(fid,'%s %f',1);
            t(n) = tmp{2};

            %skip data header line
            textscan(fid,'%s',ndheaders);

            %get data from individual file
            tmp = textscan(fid,repmat('%f ',1,ndheaders));
            fclose(fid);
            
            %grab relevant data for this flow variable
            data(:,:,:,n) = permute(reshape(tmp{q},np),order);
            
            %update user
            text = sprintf('Variable %s: %d/%d', headers{1}{q}, n, Nf);
            cpb.setValue(n/Nf*100);  	% update progress value
            cpb.setText(text);  % update user text
            
%             %extra computations for polar vectors
%             if q == 1
%                 v_tmp = permute(reshape(tmp{header_v},np),order);
%                 w_tmp = permute(reshape(tmp{header_w},np),order);
%                 
% %                 u_theta(:,:,:,n) = cos(grid.theta).*w_tmp + sin(grid.theta).*v_tmp;
% %                 u_r(:,:,:,n) = -sin(grid.theta).*w_tmp + cos(grid.theta).*v_tmp;
%                 u_r(:,:,:,n) = sin(grid.theta).*w_tmp + cos(grid.theta).*v_tmp;
%                 u_theta(:,:,:,n) = -w_tmp.*cos(grid.theta) + v_tmp.*sin(grid.theta);                
%             end
        end
        save([out,filesep,headers{1}{q},'.mat'],'data','t','grid','flags','trigger'); %need to add in actuator signal, once we get it
    end
    fprintf('\n'); %for formatting
    
%     clear data;
%     data = u_r;
%     save([out,filesep,'u_r.mat'],'data','t','grid','flags','trigger');
%     clear data;
%     data = u_theta;
%     save([out,filesep,'u_theta.mat'],'data','t','grid','flags','trigger');
end