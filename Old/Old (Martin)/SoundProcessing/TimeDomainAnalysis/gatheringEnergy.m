clear variables

mode = 'm(03)';
md = '3';

% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090921\NoiseEvents-1.5prms';
path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\20100417-CaseyData\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090917-3\NoiseEvents-1.5prms_old';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090918\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090919\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090811\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090812\NoiseEvents-1.5prms';

l = getFList(path,'mat',2,'Base');

dtUnique = (1:100)/200000;
pss_dtUnique_BL90 = zeros(length(l),length(dtUnique));
pss_dtUnique_BL30 = pss_dtUnique_BL90;

M = 1;
for n = 1:length(l)
	d = load([path '\' l{n}],'out','D','U');
	U = d.U;
	D = d.D;
	
	q = round(d.out(1).dtUnique*200000);
	if max(q) > M
		M = max(q);
	end
	pss_dtUnique_BL90(n,q) = d.out(1).pss_dtUnique;
	q = round(d.out(9).dtUnique*200000);
	if max(q) > M
		M = max(q);
	end
	pss_dtUnique_BL30(n,q) = d.out(9).pss_dtUnique;
	
end
pss_dtUnique_BL30 = mean(pss_dtUnique_BL30);
pss_dtUnique_BL90 = mean(pss_dtUnique_BL90);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = getFList(path,'mat',2,mode);

dtUnique = (1:100)/200000;
pss_dtUnique_90 = zeros(length(l),length(dtUnique));
pss_dtUnique_30 = pss_dtUnique_90;

for n = 1:length(l)
	qb = strfind(l{n},'_F');
	qe = strfind(l{n}(qb+2:end),'_');
	F(n) = str2num(l{n}(qb+2:qb+qe));
	
	d = load([path '\' l{n}],'out','D','U');
	
	q = round(d.out(1).dtUnique*200000);
	if max(q) > M
		M = max(q);
	end
	pss_dtUnique_90(n,q) = d.out(1).pss_dtUnique;
	q = round(d.out(9).dtUnique*200000);
	if max(q) > M
		M = max(q);
	end
	pss_dtUnique_30(n,q) = d.out(9).pss_dtUnique;
	
end
dtUnique = dtUnique(1:M);
[F,IX] = sort(F);
pss_dtUnique_30 = [pss_dtUnique_BL30(1:M); pss_dtUnique_30(IX,1:M)];
pss_dtUnique_90 = [pss_dtUnique_BL90(1:M); pss_dtUnique_90(IX,1:M)];

F = [0 F]*D/U*1000;
figure;
pcolor(dtUnique*U/D,F,pss_dtUnique_90); shading flat;
xlabel('\tau_j')
ylabel('St_{DF}')
title(['(Event Energy)/(Total Signal Energy) - 90^o - m = ' md])
xlim([min(dtUnique*U/D) 1+eps])
caxis([0 0.15])
colorbar
saveas(gcf,['EnergyDistribution_m' md '_90.fig'])
saveFigure_v2(gcf,['EnergyDistribution_m' md '_90'],300)
close

figure;
pcolor(dtUnique*U/D,F,pss_dtUnique_30); shading flat;
xlabel('\tau_j')
ylabel('St_{DF}')
title(['(Event Energy)/(Total Signal Energy) - 30^o - m = ' md])
xlim([min(dtUnique*U/D) 3])
caxis([0 0.07])
colorbar
saveas(gcf,['EnergyDistribution_m' md '_30.fig'])
saveFigure_v2(gcf,['EnergyDistribution_m' md '_30'],300)
close

clear IX M d l n path pss_dtUnique_BL30 pss_dtUnique_BL90 q qb qe


	
	