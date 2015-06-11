% clear variables

%Basic property looked at here is that tan(theta) = width/height. If
%tan(theta) is constant vs. forcing, then the peak is only changing height
%and the shape change is an artifact.

% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090921\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\20100417-CaseyData\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090917-3\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090918\NoiseEvents-1.5prms';
path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090919\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090811\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090812\NoiseEvents-1.5prms';

mode = 'm(03)';
md = '3';

l = getFList(path,'mat',2,'Base');

peak = [];
width = peak;
F = [];

for n = 1:length(l)
	d = load([path '\' l{n}],'out','D','U');
	U = d.U;
	D = d.D;
	BLprms(n) = d.out(9).prms;
	
	peak = [peak; abs(d.out(9).peak)];
	width = [width; d.out(9).dt];
	
end
BLprms = mean(BLprms);
BLpeak = [mean(peak(:)) std(peak(:))];
BLwidth = [mean(width(:)) std(width(:))];
BLTanTH = [mean(width(:)./peak(:)) std(width(:)./peak(:))];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = getFList(path,'mat',2,mode);

peak = zeros(length(l),2);
width = peak;
TanTH = peak;

for n = 1:length(l)
	qb = strfind(l{n},'_F');
	qe = strfind(l{n}(qb+2:end),'_');
	F(n) = str2num(l{n}(qb+2:qb+qe));
	
	d = load([path '\' l{n}],'out');
	
	ptmp = abs(d.out(9).peak);
	wtmp = d.out(9).dt;
	
	peak(n,:) = [mean(ptmp(:)) std(ptmp(:))];
	width(n,:) = [mean(wtmp(:)) std(wtmp(:))];
	TanTH(n,:) = [mean(wtmp(:)./ptmp(:)) std(wtmp(:)./ptmp(:))];
end
[F,IX] = sort(F);
width = [BLwidth; width(IX,:)];
peak = [BLpeak; peak(IX,:)];
TanTH = [BLTanTH; TanTH(IX,:)];

F = [0 F]'*D/U*1000;



clear IX M d l n path q qb qe BLTanTH BLpeak BLwidth ptmp wtmp mode

pc = polyfit(F,width(:,1),4);
wfit = polyval(pc,F);

[m,I] = max(peak(:,1));
pc = polyfit(F(1:I),peak(1:I,1),1);
pfit = polyval(pc,F);

q = find(m0.peak > m0.peak(1)*1.5,1,'last');
if ~isempty(q)
	pc = polyfit(F(I:q),peak(I:q,1),1);
	pfit(I:q) = polyval(pc,F(I:q));
	pc = polyfit(F(q+1:end),peak(q+1:end,1),4);
	pfit(q+1:end) = polyval(pc,F(q+1:end));
else
	pc = polyfit(F(I:end),peak(I:end,1),4);
	pfit(I:end) = polyval(pc,F(I:end));
end
clear m I q

eval(['m' md '.U = U;']);
eval(['m' md '.D = D;']);
eval(['m' md '.F = F;']);
eval(['m' md '.peak = peak;']);
eval(['m' md '.width = width;']);
eval(['m' md '.TanTH = TanTH;']);
eval(['m' md '.BLprms = BLprms;']);
eval(['m' md '.pfit = pfit;']);
eval(['m' md '.wfit = wfit;']);

clear md D U F TanTH pc peak pfit wfit width BLprms


% figure
% plot(m0.F,m0.peak(:,1)/m0.BLprms,'ob','MarkerFaceColor','b')
% hold on
% plot(m1.F,m1.peak(:,1)/m1.BLprms,'sr','MarkerFaceColor','r')
% plot(m3.F,m3.peak(:,1)/m3.BLprms,'d','Color',[102 204 0]/255,'MarkerFaceColor',[102 204 0]/255)
% plot(m0.F,m0.pfit/m0.BLprms,'-b','LineWidth',1.5)
% plot(m1.F,m1.pfit/m1.BLprms,'-r','LineWidth',1.5)
% plot(m3.F,m3.pfit/m3.BLprms,'Color',[102 204 0]/255,'LineWidth',1.5)
% grid on
% xlabel('St_{DF}')
% ylabel('Average Event Amplitude (p/prms_{Baseline})')
% title('Average Event Amplitude: M = 0.9, TTR = 2.0')
% legend('m = 0','m = 1','m = 3')
% 
% 
% figure
% plot(m0.F,m0.width(:,1)*1e6,'ob','MarkerFaceColor','b')
% hold on
% plot(m1.F,m1.width(:,1)*1e6,'sr','MarkerFaceColor','r')
% plot(m3.F,m3.width(:,1)*1e6,'d','Color',[102 204 0]/255,'MarkerFaceColor',[102 204 0]/255)
% plot(m0.F,m0.wfit*1e6,'-b','LineWidth',1.5)
% plot(m1.F,m1.wfit*1e6,'-r','LineWidth',1.5)
% plot(m3.F,m3.wfit*1e6,'Color',[102 204 0]/255,'LineWidth',1.5)
% grid on
% xlabel('St_{DF}')
% ylabel('Average Event Width ( \mus)')
% title('Average Event Width: M = 0.9, TTR = 2.5')
% legend('m = 0','m = 1','m = 3')
% 
% figure
% plot(m0.F,m0.wfit*1e6./m0.pfit*m0.BLprms,'-.b','LineWidth',1.5)
% hold on
% plot(m1.F,m1.wfit*1e6./m1.pfit*m1.BLprms,'--r','LineWidth',1.5)
% plot(m3.F,m3.wfit*1e6./m3.pfit*m3.BLprms,'Color',[102 204 0]/255,'LineWidth',1.5)
% grid on
% xlabel('St_{DF}')
% legend('m = 0','m = 1','m = 3')
% ylabel('Shape Ratio: width/height')
% title('Average Event Shape Ratio: M = 0.9, TTR = 1.0')
