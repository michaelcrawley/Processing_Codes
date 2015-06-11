clear;

%Constants
phys.Uj = 285.99; 
phys.D = 0.0254;
dt = 0.001*phys.D/phys.Uj; %will need to verify that files aren't missing, and then calculate the FS
phys.M = 0.9;
phys.Temp = 251.31; %static jet temperature
matversion = 1.21;

%The snapshot files will need to be grabbed as an excitation group, and
%array
[flist,src] = getfiles('*.mat',[]);
fname = input('What is the desired filename? ','s');
out_dir = uigetdir('Specify output directory',pwd);

%grab iteration numbers from filenames, reorder flist
itr = zeros(size(flist));
for n = 1:length(flist)
    tmp = regexp(flist{n},'\d+','match');
    itr(n) = str2double(tmp{1});
end
[itr,I] = sort(itr);
d_itr = unique(diff(itr));
if length(d_itr) > 1
    error('Missing files!!!');
end
flist = flist(I);

%Read in grid data
gridfile = getfiles('*.grd',src); 
griddata = dlmread([src filesep gridfile{1}]); %output is in i,j,x,y,z columns and is normalized by D
x = griddata(:,1);
y = sqrt(griddata(:,2).^2 + griddata(:,3).^2);

%Calculate additional constants
phys.To = phys.Temp*(1 + .2*phys.M^2);
phys.a = sqrt(1.4*287*phys.To);
phys.FS = 1/(dt*d_itr);

%Grab data
p = zeros(length(flist),length(x));
for n = 1:length(flist)
    tmp = load([src filesep flist{n}]);
    p(n,:) = tmp.temp;
end

%subtract to get fluctuations
pmn = mean(p,1);
p = p-repmat(pmn,length(flist),1);

%Interpolate onto a regular axial grid
mndx = mean(diff(x));
chk = [0; abs(mndx-diff(x)) > 1e-4];
nf.lblocks.p = p(:,~chk);
phys.x = x(~chk);
phys.y = y(~chk);

save([out_dir filesep fname],'nf','phys','matversion');