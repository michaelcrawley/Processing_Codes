% clear variables

W = 1.23; phD = [1/4 1/2 3/4 1/4 1/2 1/4];   %expected wavelength (x/D)

St = 0.52;
Uj = 477;

x0 = 0;

clear d;
IM = vm11_St052.CondAvg.Q; dx = mean(diff(m0_St013.x(:,1)));
yl = find(m0_St013.y(1,:) < 0.75,1,'first');
yl(2) = find(m0_St013.y(1,:) < -0.75,1,'first');
IM = permute(IM(:,yl(1):yl(2),:),[2 1 3]); clear yl;

D = 0.0254; %m - jet diameter
ainf = sqrt(1.4*287.05*(273.15+25));

X = [];
V = [];

[ds,corr,c,vn,vp] = PIVCorrelate(IM(:,:,1),IM(:,:,2),W/dx,phD(1));
X(:,end+1) = (vp(:,1)-x0)*dx;  %x/D streamwise distance
V(:,end+1) = vp(:,2)*dx*St*Uj*4;   %m/s - Convert pixel displacement into velocity
X(:,end+1) = (vn(:,1)-x0)*dx;
V(:,end+1) = -vn(:,2)*dx*St*Uj*4/3;
disp('2')
[ds,corr,c,vn,vp] = PIVCorrelate(IM(:,:,1),IM(:,:,3),W/dx,phD(2));
X(:,end+1) = (vp(:,1)-x0)*dx;
V(:,end+1) = vp(:,2)*dx*St*Uj*2;
X(:,end+1) = (vn(:,1)-x0)*dx;
V(:,end+1) = -vn(:,2)*dx*St*Uj*2;
disp('3')
[ds,corr,c,vn,vp] = PIVCorrelate(IM(:,:,1),IM(:,:,4),W/dx,phD(3));
X(:,end+1) = (vp(:,1)-x0)*dx;
V(:,end+1) = vp(:,2)*dx*St*Uj*4/3;
X(:,end+1) = (vn(:,1)-x0)*dx;
V(:,end+1) = -vn(:,2)*dx*St*Uj*4;

disp('4')
[ds,corr,c,vn,vp] = PIVCorrelate(IM(:,:,2),IM(:,:,3),W/dx,phD(4));
X(:,end+1) = (vp(:,1)-x0)*dx;
V(:,end+1) = vp(:,2)*dx*St*Uj*4;
X(:,end+1) = (vn(:,1)-x0)*dx;
V(:,end+1) = -vn(:,2)*dx*St*Uj*4/3;
disp('5')
[ds,corr,c,vn,vp] = PIVCorrelate(IM(:,:,2),IM(:,:,4),W/dx,phD(5));
X(:,end+1) = (vp(:,1)-x0)*dx;
V(:,end+1) = vp(:,2)*dx*St*Uj*2;
X(:,end+1) = (vn(:,1)-x0)*dx;
V(:,end+1) = -vn(:,2)*dx*St*Uj*2;

disp('6')
[ds,corr,c,vn,vp] = PIVCorrelate(IM(:,:,3),IM(:,:,4),W/dx,phD(6));
X(:,end+1) = (vp(:,1)-x0)*dx;
V(:,end+1) = vp(:,2)*dx*St*Uj*4;
X(:,end+1) = (vn(:,1)-x0)*dx;
V(:,end+1) = -vn(:,2)*dx*St*Uj*4/3;
clear ds corr c vn vp

figure; hold on;
CC = {'or','or','og','og','ob','ob','sr','sr','sg','sg','dr','dr'};
% CC = {'or','or','og','og','sr','sr'};
for n = 1:size(V,2)
    plot(X(:,n),V(:,n)/ainf,CC{n});
end

x = unique(X(:));
v = zeros(size(x)); verr = v;
q = V(6:end-6,:); qm = mean(q(:)); qs = std(q(:));
qw = [phD; 1-phD]; qw = 1./qw(:)'; qw = repmat(qw,[size(V,1) 1]);
for n = 1:length(x)
    I = logical((X==x(n)).*(V~=0).*((V >= qm-3*qs) & (V <= qm+3*qs)));
    q = V(I); 
    qwt = qw(I);
    if ~isempty(q)
        v(n) = sum(q(:).*qwt(:))/sum(qwt(:));
        verr(n) = sqrt(sum(qwt(:).*(q(:)-v(n)).^2)/sum(qwt(:)));
    end
end
clear n I q qw qwt CC
errorbar(x,v/ainf,verr/ainf,'-k','LineWidth',1.5)
axis([1 12 0 1.4+eps])
xlabel('x/D')
ylabel('U_c/a_\infty')
grid on


%% EDIT ME %%%%%%%%%%%%%%%%%%%
title('M = 1.3, T_o/T_a = 1.5, m = \pm1, St_{DF} = 0.52')
% save Uc-M13_TTR10_m0_F11299.mat
saveas(gcf,'Uc-M13_TTR15_vm11_St052.fig')
saveFigure_v2(gcf,'Uc-M13_TTR15_vm11_St052',600)
