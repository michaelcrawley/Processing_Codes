clear variables

% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090921\NoiseEvents-1.5prms';
path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\20100417-CaseyData\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090917-3\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090918\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090919\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090811\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090812\NoiseEvents-1.5prms';

mode = 'm(00)';
md = '0';

l = getFList(path,'mat',2,'Base');

t = (-499:500)/200000;
avg_BL90 = zeros(length(l),length(t));
avg_BL30 = avg_BL90;

for n = 1:length(l)
	d = load([path '\' l{n}],'out','D','U');
	U = d.U;
	D = d.D;
	
	avg_BL90(n,:) = d.out(1).avg;
	avg_BL30(n,:) = d.out(9).avg;
	
end
avg_BL30 = mean(avg_BL30);
avg_BL90 = mean(avg_BL90);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = getFList(path,'mat',2,mode);

avg_90 = zeros(length(l),length(t));
avg_30 = avg_90;

for n = 1:length(l)
	qb = strfind(l{n},'_F');
	qe = strfind(l{n}(qb+2:end),'_');
	F(n) = str2num(l{n}(qb+2:qb+qe));
	
	d = load([path '\' l{n}],'out','D','U');
	
	avg_90(n,:) = d.out(1).avg;
	avg_30(n,:) = d.out(9).avg;
	
end
[F,IX] = sort(F);
avg_30 = [avg_BL30; avg_30(IX,:)];
avg_90 = [avg_BL90; avg_90(IX,:)];

F = [0 F]*D/U*1000;
% figure;
% pcolor(t*U/D,F,avg_90); shading flat;
% xlabel('\tau_j')
% ylabel('St_{DF}')
% title(['Average Signal Shape - 90^o - m = ' md])
% xlim([-2 2])
% caxis([-0.1 0.4])
% colorbar
% saveas(gcf,['ShapeDistribution_m' md '_90.fig'])
% saveFigure_v2(gcf,['ShapeDistribution_m' md '_90'],300)
% close
% 
% figure;
% pcolor(t*U/D,F,avg_30); shading flat;
% xlabel('\tau_j')
% ylabel('St_{DF}')
% title(['Average Signal Shape - 30^o - m = ' md])
% xlim([-8 8])
% caxis([-0.5 0.6+eps])
% colorbar
% saveas(gcf,['ShapeDistribution_m' md '_30.fig'])
% saveFigure_v2(gcf,['ShapeDistribution_m' md '_30'],300)
% close


clear IX M d l n path avg_BL30 avg_BL90 q qb qe


	
	