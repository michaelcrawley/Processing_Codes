function [x,y,z,u,v,w,c] = read_vel_data(xshift,yshift,To)
%This function reads a standard exported PIV vector sets into matrices when
%given the horizontal and vertical
%shifts (in mm) to apply to the data. Note: Shifts are positive when
%shifting the zero up and to the right, and the stagnation temperature of
%the flow

global To;

dir_name = uigetdir('My Computer','Specify Directory Containing Data'); 
flist = struct2cell(dir(dir_name));
q = strfind(flist(1,:),'.dat'); %extracts list of data files - ignores all other files
keep = zeros(size(q));
for n = 1:length(keep)
    if ~isempty(q{n})
        keep(n) = 1;
    end
end
flist = flist(1,logical(keep));
[keep,ok] = listdlg('PromptString','Select files to process:',...
                'SelectionMode','multiple',...
                'ListString',flist);    %Asks the user to select the subset of data files they wish to process
if ok==0
    error('PT_UserCancel','Program Terminated Due to User Selection of Cancel')
end
flist = flist(keep);
fn=flist{1};


fid=fopen([dir_name,'\',fn]);

if fid>0
    fgets(fid);
    fgets(fid);
    
    fscanf(fid,'%c',20);
    kmax = fscanf(fid,'%u',2);
    fscanf(fid,'%c',4);
    jmax = fscanf(fid,'%u',2);
    fgets(fid);
    data = fscanf(fid,'%g %g %g %g %g %g',[6 inf]);

    xi(:,:) = reshape(data(1,:),kmax,jmax)-xshift;
    yi(:,:) = reshape(data(2,:),kmax,jmax)-yshift;
    zi(:,:) = reshape(data(3,:),kmax,jmax);
    
    ui(:,:) = reshape(data(6,:),kmax,jmax);
    vi(:,:) = reshape(data(5,:),kmax,jmax);
    wi(:,:) = reshape(data(4,:),kmax,jmax);

    for l = size(yi,2):-1:1
        if yi(1,l)<0 & yi(1,l-1) >= 0
            ipy=l-1;
            break
        end
    end

    ipx = 5;
    epx= 237;
    epy= 10;

    x(:,:)=xi(ipx:epx,ipy:-1:epy);
    y(:,:)=yi(ipx:epx,ipy:-1:epy);
    z(:,:)=zi(ipx:epx,ipy:-1:epy);

    u(:,:)=ui(ipx:epx,ipy:-1:epy);
    v(:,:)=vi(ipx:epx,ipy:-1:epy);
    w(:,:)=wi(ipx:epx,ipy:-1:epy);

    c(:,:)=sqrt(1.4*287*(To-(u(:,:).^2+v(:,:).^2+w(:,:).^2)/(2*1004.5)));
end

fclose(fid);