function [x,z,u,w] = readDAT(fname)
    %reads 2D DAT files

fid = fopen(fname,'r');
txt = fgetl(fid); txt = fgetl(fid);
fscanf(fid, '%c',7);
Imax = fscanf(fid, '%u2');
fscanf(fid, '%c', 4);
Jmax = fscanf(fid, '%u2');
fclose(fid);

% q = strfind(txt,'I=');
% qq = strfind(txt,',');
% q2 = qq(min(find(qq>q)));
% Imax = str2num(txt(q+2:q2-1));
% q = strfind(txt,'J=');
% qq = strfind(txt,',');
% q2 = qq(min(find(qq>q)));
% Jmax = str2num(txt(q+2:q2-1));

data = dlmread(fname,' ',3,0);
x = reshape(data(:,1),Imax,Jmax);
z = reshape(data(:,2),Imax,Jmax);
u = reshape(data(:,3),Imax,Jmax);
w = reshape(data(:,4),Imax,Jmax);
