function varargout = readPIVfiles_IMX_v1()
%program processes all files in a directory
%ALL FILES MUST BE SAME SIZE

D = 0.0254;
M = 1.65;
Ta = 273.15+20.5;   %ambient temperature - K
SD = 'F:\Shared Docs\PIV\20090212';
TAG = 'M165';

fdir = uigetdir(SD,'Specify Directory Containing Data');
if fdir==0 
    error('Program Terminated Due to No Directory Selection')
end

To = input('Enter Stagnation Temp (C): ')+273.15;    %measured stagnation temp - K

Tj = To/(1+1/5*M^2);    %jet exit temp - K

flist = struct2cell(dir(fdir));
flist = flist(1,:);
keep = [];
for n = 1:length(flist)
    [pathstr, name, ext]=fileparts(flist{n});
    if strcmpi(ext,'.vc7')
        keep = [keep n];
    end
end
flist = flist(keep);
fk = [];
for n = 1:length(flist)
    progress(n,1,length(flist),10);
    if n==1
        RDAT = readimx(fullfile(fdir,flist{n}));
        [x,y,U,V] = showimx(RDAT);
        x = x/1000/D;
        y = y/1000/D;
    else
        RDAT = readimx(fullfile(fdir,flist{n}));
        [xT,yT,U(:,:,n),V(:,:,n)] = showimx(RDAT);
    end
end
% U = U(:,:,fk);
% V = V(:,:,fk);
Um = mean(U,3);
Vm = mean(V,3);
% surface(x,y,sqrt(Um.^2+Vm.^2),'EdgeColor','none');
% shading interp;
% axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])

x = x(6:end,:); y = y(6:end,:); U = U(6:end,:,:); V = V(6:end,:,:); Um = Um(6:end,:); Vm = Vm(6:end,:);

qx = min(find(x(:,1)>1));
qy = max(find(y(1,:)>0));
jet = Um(1:qx,qy-1:qy+1);
Uj = mean(jet(:));  %jet exit velocity - m/s

CL = Um(:,qy)/sqrt(7/5*287.05*Tj);  %jet centerline mach number - dimensionless
TKE = 1/Uj^2*(std(U,0,3).^2 + std(V,0,3).^2);  %2D TKE - dimensionless
FWHM = findFWHM(x,y,Um);    %jet width - y/D

qx = strfind(fdir,TAG); qx = qx(1);
qx2 = strfind(fdir,'\');
qx2 = qx2(qx2 > qx); qx2 = qx2(1);

if ~exist([fdir(1:qx-1) 'Processed'],'dir')
    mkdir([fdir(1:qx-1) 'Processed']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % qq = Um(5:7,49:52);
% % % TKE = TKE*Uj^2/mean(qq(:))^2;
% % % Uj = mean(qq(:)); clear qq
% % % 
% % % CL = CL(5:end);
% % % FWHM = FWHM(5:end);
% % % TKE = TKE(5:end,:);
% % % U = U(5:end,:,:);
% % % V = V(5:end,:,:);
% % % Um = Um(5:end,:);
% % % Vm = Vm(5:end,:);
% % % x = x(5:end,:);
% % % y = y(5:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([fdir(1:qx-1) 'Processed\' fdir(qx:qx2-1) '.mat'],'x','y','U','V','Um','Vm','Uj','To','Tj','CL','TKE','FWHM')

% qx = strfind(fdir,'Raw');
% save([fdir(1:qx(1)-1) 'Processed\reorgData.mat'],'x','y','U','V','Um','Vm','Uj','To','Tj','CL','TKE','FWHM')
    
if nargout > 0
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = U;
    varargout{4} = V;
    varargout{5} = Um;
    varargout{6} = Vm;
    varargout{7} = Uj;
    varargout{8} = To;
    varargout{9} = Tj;
    varargout{10} = CL;
    varargout{11} = TKE;
    varargout{12} = FWHM;
else
    varargout{1} = 'No Output';
end