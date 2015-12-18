function v = loadpivtxt(fname)
%LOADPIVTXT  Load a vector field exported in text format
%   V = LOADPIVTXT(FILENAME) loads a vector field FILENAME saved in text
%   format (TXT or DAT) into the structure V. Accepted file formats are:
%     .TXT, exported from DaVis (LaVision)
%     .DAT, from VidPIV (Oxford Laser)
%
%   The fields of the structure V are the same as for LOADVEC, except the
%   field V.Attributes, which contains the first line of the text file
%   FILENAME.
%
%   Using LOADVEC to load vector fields in TXT mode allows for wildcards
%   and brackets (filename expansion).
%
%   Note: Using VEC/VC7 files instead of TXT files is strongly encouraged.
%   If you prefer platform-independant file formats, you may convert your
%   original DaVis files into Mat-files using VEC2MAT.
%
%   Example:
%      v=loadpivtxt('B00001.txt');
%
%   See also LOADVEC, VEC2MAT.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.32,  Date: 2012/07/11
%   This function is part of the PIVMat Toolbox


% History:
% 2002/30/01: v1.00, first version.
% 2003/02/17: v1.01.
% 2005/02/23: v1.02, cosmetic changes.
% 2005/09/01: v1.03, argument check added.
% 2005/10/06: v1.04, file validity checked.
% 2005/10/11: v1.10, loadpivtxt is combination of readpivtxt v1.04 and
%                    dav2mat v1.05 (2005/10/06).
% 2005/10/17: v1.11, bug fixed for vectors x and y.
% 2005/10/21: v1.20, bug Y-axis up/down-ward fixed. Reads the units in the
%                    commentstring (1st line)
% 2006/03/03: v1.21, history is now a cell array of strings.
% 2006/09/12: v1.22, new field 'source'
% 2010/09/27: v1.30, also work when the comment line is missing
% 2012/06/22: v1.31, help text changed: accept VidPIV files
% 2012/07/11: v1.32, bug fixed regarding the scan of the comment line


error(nargchk(1,1,nargin));

fid=fopen(fname);
if fid<0,
    error(['Can''t open file ' fname]);
else
    commentstring = fscanf(fid,'%s',12);
    if ~strfind(commentstring,'#')    %new v1.30 : if the comment is missing
        commentstring = '#missingcomment "" "mm" "1" "mm" "velocity" "m/s"';
        fclose(fid);
        fid=fopen(fname);
    end
    [data,count]= fscanf(fid,'%f',[4,inf]);
    fclose(fid);
    if ~count %new 1.04
        error(['''' fname ''' is not a valid .TXT or .DAT file']);
    end
    
    nx=length(find(data(2,:)==data(2,1)));
    ny=size(data,2)/nx;
    
    v.x=zeros(nx,ny);
    v.y=zeros(nx,ny);
    v.vx=zeros(nx,ny);
    v.vy=zeros(nx,ny);
    
    for j=1:ny,
        for i=1:nx,
            v.x(i,j)=data(1,(j-1)*nx+i);
            v.y(i,j)=data(2,(j-1)*nx+i);
            v.vx(i,j)=data(3,(j-1)*nx+i);
            v.vy(i,j)=data(4,(j-1)*nx+i);
        end
    end
    
    v.x=v.x(:,1)';
    v.y=v.y(1,:);
    
    % new v1.20
    p=findstr(commentstring,'"');
    v.unitx=strtok(commentstring((p(3)+1):end),'"');
    v.unity=strtok(commentstring((p(7)+1):end),'"');
    v.unitvx=strtok(commentstring((p(11)+1):end),'"');
    v.unitvy=strtok(commentstring((p(11)+1):end),'"');
    
    v.namevx='u_x';
    v.namevy='u_y';
    v.namex='x';
    v.namey='y';
    
    % new v1.20:
    if v.y(2)>v.y(1),
        % natural orientation for DaVis (mode axis IJ (matrix)):
        v.ysign='Y axis downward';
    else
        % anti-natural orientation for DaVis (mode axis XY):
        v.ysign='Y axis upward';
        % inverts up/down:
        v.vx =  v.vx(:, end:-1:1);
        v.vy =  v.vy(:, end:-1:1); % No nead to invert the vy component in txt mode
        v.y = v.y(end:-1:1); % the vector y must have increasing components
    end
    
    v.setname=getsetname;
    v.name=fname;
    v.Attributes=commentstring;
    v.history={'loadpivtxt'};
    if strfind(lower(fname),'.txt')
        v.source='TXT file';
    elseif strfind(lower(fname),'.dat')
        v.source='DAT file';
    else
        v.source='Text file';
    end
end
