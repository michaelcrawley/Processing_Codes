clear variables

% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090921\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\20100417-CaseyData\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090917-3\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090918\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090919\NoiseEvents-1.5prms';
path = 'F:\Research\Mach1.3\Acoustics\Array\20090811\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090812\NoiseEvents-1.5prms';

mode = 'm(03)';
md = '3';

l = getFList(path,'mat',2,'Base');

dtHist_BL90 = zeros(length(l),10);
dtHist_BL30 = dtHist_BL90;

for n = 1:length(l)
	d = load([path '\' l{n}],'dtHist','D','U');
	U = d.U;
	D = d.D;
	
	tmp = d.dtHist(:,1)'; dtHist_BL90(n,1:length(tmp)) = tmp;
	tmp = d.dtHist(:,9)'; dtHist_BL30(n,1:length(tmp)) = tmp;
	
end
dtHist_BL30 = mean(dtHist_BL30);
dtHist_BL90 = mean(dtHist_BL90);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = getFList(path,'mat',2,mode);

dtHist_90 = zeros(length(l),10);
dtHist_30 = dtHist_90;

for n = 1:length(l)
	qb = strfind(l{n},'_F');
	qe = strfind(l{n}(qb+2:end),'_');
	F(n) = str2num(l{n}(qb+2:qb+qe));
	
	d = load([path '\' l{n}],'dtHist','D','U');
	
	tmp = d.dtHist(:,1)'; dtHist_90(n,1:length(tmp)) = tmp;
	tmp = d.dtHist(:,9)'; dtHist_30(n,1:length(tmp)) = tmp;
	
end
[F,IX] = sort(F);
M = min([length(dtHist_BL30) size(dtHist_30,2)]);
dtHist_30 = [dtHist_BL30(1:M); dtHist_30(IX,1:M)];
dtHist_90 = [dtHist_BL90(1:M); dtHist_90(IX,1:M)];

F = [0 F]*D/U*1000;
% figure;
% pcolor((1:size(dtHist_90,2))/200000*U/D,F,log10(dtHist_90)); shading flat;
% % mesh((1:size(dtHist_90,2))/200000*U/D,F,dtHist_90);
% % set(gca,'ZScale','log')
% % view(15,30)
% xlabel('\tau_j')
% ylabel('St_{DF}')
% title(['Peak Width Histogram - 90^o - m = ' md])
% xlim([0 4])
% caxis([-5 -0.5])
% colorbar
% saveas(gcf,['WidthDistribution_m' md '_90.fig'])
% saveFigure_v2(gcf,['WidthDistribution_m' md '_90'],300)
% close

% figure;
% pcolor((1:size(dtHist_90,2))/200000*U/D,F,log10(dtHist_30)); shading flat;
% % mesh((1:size(dtHist_90,2))/200000*U/D,F,dtHist_30);
% % set(gca,'ZScale','log')
% % view(15,30)
% xlabel('\tau_j')
% ylabel('St_{DF}')
% title(['Peak Width Histogram - 30^o - m = ' md])
% xlim([0 4])
% caxis([-5 -0.5])
% colorbar
% saveas(gcf,['WidthDistribution_m' md '_30.fig'])
% saveFigure_v2(gcf,['WidthDistribution_m' md '_30'],300)
% close

% figure; hold on;
% CM = colormap('jet'); CM(1,:) = [0 0 0];
% for n = 1:29
% 	plot((1:size(dtHist_90,2))/200000*U/D,dtHist_30(n,:),'Color',CM(floor(64/29*n),:))
% end
% set(gca,'YScale','log')
% xlim([0 3])

clear IX M d l n path dtHist_BL30 dtHist_BL90 q qb qe tmp


	
	