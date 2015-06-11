clear variables

mode = 'm(00)';
md = '0';
CH = 9;
ANG = '90';

% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090921\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\20100417-CaseyData\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090917-3\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090918\NoiseEvents-1.5prms';
path = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090919\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090811\NoiseEvents-1.5prms';
% path = 'F:\Research\Mach1.3\Acoustics\Array\20090812\NoiseEvents-1.5prms';

l = getFList(path,'mat',2,'Base');

dt0x = (1:5000);	%initial time steps for histogram

M = 0;
dt0All = zeros(10000,length(l));
NB = length(l);	%number of baselines
for n = 1:NB
	d = load([path '\' l{n}],'out','D','U');
	U = d.U;
	D = d.D;
	
	t0 = d.out(CH).t0;
	pk = d.out(CH).peak;
% 	pkl = pk > 0;	%Select positive peaks only
	pkl = true(size(pk));	%Select all peaks

	t0p = t0(pkl);

	N = 1; %order of the most distant neighbor
	dt0 = zeros(length(t0p)-2*N,N);
	for m = 1:N		%Calculates differences and converts into an integer number of samples
		dt0(:,m) = round((t0p(N+1:end-N) - t0p(N+1-m:end-N-m))*200000);
	end
	if max(dt0) > M
		M = max(dt0(:));
	end
	
	Nevents(n) = length(dt0(:,N));
	dt0All(1:length(dt0),n) = dt0(:,N);	
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = getFList(path,'mat',2,mode);

for n = 1:length(l)
	qb = strfind(l{n},'_F');
	qe = strfind(l{n}(qb+2:end),'_');
	F(n) = str2num(l{n}(qb+2:qb+qe));
	
	d = load([path '\' l{n}],'out');
	
	t0 = d.out(CH).t0;
	pk = d.out(CH).peak;
% 	pkl = pk > 0;	%Select positive peaks only
	pkl = true(size(pk));	%Select all peaks

	t0p = t0(pkl);

	N = 1; %order of the most distant neighbor
	dt0 = zeros(length(t0p)-2*N,N);
	for m = 1:N		%Calculates differences and converts into an integer number of samples
		dt0(:,m) = round((t0p(N+1:end-N) - t0p(N+1-m:end-N-m))*200000);
	end
	if max(dt0) > M
		M = max(dt0(:));
	end
	
	Nevents(n+NB) = length(dt0(:,N));
	dt0All(1:length(dt0),n+NB) = dt0(:,N);
end

dt0x = dt0x(1:M);	%eliminate all time steps larger than needed

dt0Hist = hist(dt0All,dt0x);% dt0Hist = dt0Hist./repmat(Nevents,size(dt0Hist,1),1);
dt0x = dt0x(2:end); dt0Hist = dt0Hist(2:end,:);	%Discards zero events since they are artifacts of empty cells in the matrix
dt0x = dt0x/200000*U/D;	%Converts distribution points into units of inverse Strouhal number

dt0Hist(:,NB) = mean(dt0Hist(:,1:NB),2); dt0Hist = dt0Hist(:,NB:end);	%averages the histograms of the baselines

F = [0 F]*D/U*1000;
[F,IX] = sort(F); dt0Hist = dt0Hist(:,IX)';

figure;
pcolor(dt0x,F,dt0Hist); shading flat;
xlabel('\tau_j')
ylabel('St_{DF}')
title(['(# of Intervals of Given Duration)/(Total # of intervals) - ' ANG '^o - m = ' md])
xlim([min(dt0x) 10])
caxis([0 0.04])
colorbar
saveas(gcf,['IntervalDistribution_m' md '_' ANG '-pospeaks.fig'])
saveFigure_v2(gcf,['IntervalDistribution_m' md '_' ANG '-pospeaks'],300)
close

% load temp1.mat
% tmp(CH).x = dt0x; tmp(CH).h = dt0Hist;
% save('temp1.mat','tmp')

	
%%%%%%%%%%% PROTOTYPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % for CH = 1:10
% % % 	t0 = out(CH).t0;
% % % 	pk = out(CH).peak;
% % % 	pkl = pk > 0;	%Select positive peaks only
% % % 	% pkl = true(size(pk));	%Select all peaks
% % % 
% % % 	t0p = t0(pkl); Nevents(CH) = length(t0p);
% % % 
% % % 	N = 5; %order of the most distant neighbor
% % % 	dt0 = zeros(length(t0p)-2*N,N);
% % % 	for m = 1:N		%Calculates differences and converts into an integer number of samples
% % % 		dt0(:,m) = round((t0p(N+1:end-N) - t0p(N+1-m:end-N-m))*200000);
% % % 	end
% % % % 	dt0x = unique(dt0);
% % % % 
% % % % 	dt0Hist = hist(dt0,dt0x); dt0Hist = dt0Hist/length(t0p);
% % % % 	dt0x = dt0x/200000*U/D;	%Converts distribution points into units of inverse Strouhal number
% % % % 
% % % % 	figure;	%Plot in units of inverse Strouhal number with normalization for spacing size
% % % % 	plot(dt0x,dt0Hist(:,1))
% % % % 	hold on
% % % % 	plot(dt0x/2,dt0Hist(:,2)*2,'g')
% % % % 	plot(dt0x/3,dt0Hist(:,3)*3,'r')
% % % % 	plot(dt0x/4,dt0Hist(:,4)*4,'k')
% % % % 	plot(dt0x/5,dt0Hist(:,5)*5,'m')
% % % % 	xlim([0 10])
% % % 
% % % 
% % % 	dt0All(1:length(dt0(:,1)),CH) = dt0(:,1);
% % % end
% % % dt0x = unique(dt0All);
% % % 
% % % dt0Hist = hist(dt0All,dt0x); dt0Hist = dt0Hist./repmat(Nevents,size(dt0Hist,1),1);
% % % dt0x = dt0x(2:end); dt0Hist = dt0Hist(2:end,:);	%Discards zero events since they are artifacts of empty cells in the matrix
% % % dt0x = dt0x/200000*U/D;	%Converts distribution points into units of inverse Strouhal number
% % % figure; CM = colormap('jet'); hold on;
% % % for n = 1:10
% % % 	plot(dt0x,dt0Hist(:,n),'Color',CM(floor(6.4*n),:))
% % % end
% % % xlim([0 10])
% % % legend('90','80','70','60','50','45','40','35','30','25')
