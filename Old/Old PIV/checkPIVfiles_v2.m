function [x,y,U,V,Um,Vm,CL,TKE,Uj,To,Tj] = checkPIVfiles_v2()
%program processes all files in a directory
%ALL FILES MUST BE SAME SIZE

D = 0.0254;
M = 1.3;
% Po = 9.9;   %flow stagnation pressure - psig
% Paux = 70;  %seeder auxiliary pressure - psig
Ta = 273.15+20;   %ambient temperature - K

fdir = uigetdir('Specify Directory Containing Data',cd);
if fdir==0 
    error('Program Terminated Due to No Directory Selection')
end

q = strfind(fdir,'C_');
q2 = strfind(fdir,'_');
qq = max(find(q2<q));
Tom = str2num(fdir(q2(qq)+1:q-1))+273.15;    %measured stagnation temp - K

    %corrects stagnation temperature for temperature drop introduced by
    %solid particle seeder -NOT SURE THIS WORKS CORRECTLY
% ms = (6894.76*(Paux-Po)+1.01e5)/(1.577*287.05*Ta)*pi*(0.0254/16)^2*sqrt(1.4*287.05*Ta/1.2);   %seeder mass flow rate
% mt = (6894.76*Po+1.01e5)/(1.456*287.05*Tom)*pi*(0.0254/2)^2*sqrt(1.4*287.05*Tom/1.162);  %total mass flow rate
% To = Tom*(1-ms/mt)+Ta*ms/mt; %average temperature - K
To = Tom;

Tj = To/(1+1/5*M^2);    %jet exit temp - K

flist = struct2cell(dir(fdir));
flist = flist(1,:);
keep = [];
for n = 1:length(flist)
    [pathstr, name, ext]=fileparts(flist{n});
    if strcmp(ext,'.dat')
        keep = [keep n];
    end
end
flist = flist(keep);
fk = [];
for n = 1:length(flist)
    progress(n,1,length(flist),10);
    if n==1
        [x,y,U,V] = readDAT(fullfile(fdir,flist{n}));
        x = x/1000/D;
        y = y/1000/D;
        
        surface(x,y,sqrt(U.^2+V.^2),'EdgeColor','none');
        shading interp;
        axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])
        title(flist{n})
        pause(0.7)
        close
%         hold on
%         quiver3(x,y,ones(size(x))*max(sqrt(U(:).^2+V(:).^2)),U,V,zeros(size(x)),'w')
%         hold off
%         axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])
%         colorbar
%         kf = input('Keep File: ','s');
%         close
%         if ~strcmp(kf,'')
%             delete(fullfile(fdir,flist{n}))
%         else
%             fk = [fk n];
%         end
    else
        [x,y,U(:,:,n),V(:,:,n)] = readDAT(fullfile(fdir,flist{n}));
        x = x/1000/D;
        y = y/1000/D;
        
%         surface(x,y,sqrt(u(:,:,n).^2+V(:,:,n).^2),'EdgeColor','none');
%         shading interp;
%         axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])
%         title(flist{n})
%         pause(0.7)
%         close

%         hold on
%         quiver3(x,y,ones(size(x))*max(max(sqrt(U(:,:,n).^2+V(:,:,n).^2))),U(:,:,n),V(:,:,n),zeros(size(x)),'w')
%         hold off
%         axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])
%         colorbar
%         kf = input('Keep File: ','s');
%         close
%         if ~strcmp(kf,'')
%             delete(fullfile(fdir,flist{n}))
%         else
%             fk = [fk n];
%         end
    end
end
% U = U(:,:,fk);
% V = V(:,:,fk);
Um = mean(U,3);
Vm = mean(V,3);
surface(x,y,sqrt(Um.^2+Vm.^2),'EdgeColor','none');
shading interp;
axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])

qx = min(find(x(:,1)>1));
qy = max(find(y(1,:)>0));
jet = U(1:qx,qy-1:qy+1);
Uj = mean(jet(:));  %jet exit velocity - m/s

CL = Um(:,qy)/sqrt(7/5*287.05*Tj);  %jet centerline velocity - dimensionless
TKE = 1/Uj^2*(std(U,0,3).^2 + std(V,0,3).^2);  %2D TKE - dimensionless

    