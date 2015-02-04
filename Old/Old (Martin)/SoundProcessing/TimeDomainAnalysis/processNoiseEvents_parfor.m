%%%% Processes acoustic data for time-domain analysis and wavelet
%%%% reconstruction. Uses parallel processing.

% pth = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090921\';
pth = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\20100417-CaseyData\';
% pth = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090917-3\';
% pth = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090918\';
% pth = 'F:\Research\Mach0.9\ForcingNoiseData\Array Measurements\Recollecting 200909\20090919\';
% pth = 'F:\Research\Mach1.3\Acoustics\Array\20090811\';
% pth = 'F:\Research\Mach1.3\Acoustics\Array\20090812\';

	%Threshold for peak identification is 1.5prms of the baseline for all
	%signals.
% FL = getFList(pth,'NOS',2,'Base');
% p = [];
% for n = 1:length(FL)
% 	d = load([pth FL{n}]);
% 	p(end+1,:) = std(d);
% end
% thresh = mean(p)*1.5; clear p;

FL = getFList(pth,'NOS',2,'Base');	%List of files to process.
FL = cat(2,FL,getFList(pth,'NOS',2,'m(00)'));
% FL = cat(2,FL,getFList(pth,'NOS',2,'m(01)(01)'));
% FL = cat(2,FL,getFList(pth,'NOS',2,'m(03)'));

	%Jet exit diameter and velocity.
D = 0.0254;
U = Vj(0.9,1.0,300); U = U.Vj;

Nevents = zeros(length(FL),10);	%Stores the number of events that occurred in each data stream
matlabpool open 4;	%initiates parallel processing
for m = 1:length(FL)
	F = FL{m};

	disp(['..' F])
	d = dlmread([pth F],'\t');
	thresh = std(d)*1.5;	%Threshold for peak identification is 1.5prms of the signal being analyzed.

	if ~exist([pth 'NoiseEvents-1.5prms'],'dir')
		mkdir(pth,'NoiseEvents-1.5prms\dtHist')
		mkdir(pth,'NoiseEvents-1.5prms\Energy')
		mkdir(pth,'NoiseEvents-1.5prms\Shape')
	end

	disp('...Processing Events...')
	[signals,out] = AccSignal_v2_parfor(d,thresh);	%performs signal analysis
	for n = 1:10
		Nevents(m,n) = length(out(n).dt);
	end

	%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure; hold on;
	CM = colormap('jet');
	for n = 1:10
		plot(out(n).dtUnique*U/D,out(n).pss_dtUnique,'-o','Color',CM(n*6,:))
	end
	legend('90','80','70','60','50','45','40','35','30','25')
	axis([0 4 0 0.25+eps])
	grid on
	xlabel('\tau_j')
	ylabel('(Event Energy)/(Total Signal Energy)')
	title('Normalized Energy of Events of Various Widths')
	saveFigure_v2(gcf,[pth 'NoiseEvents-1.5prms\Energy\' F(1:end-4) '_Energy'],300)
	close

	figure; hold on;
	for n = 10:-1:1
		plot((-499:500)'/200000*U/D,out(n).avg/out(n).prms,'-','Color',CM(n*6,:))
	end
	legend('25','30','35','40','45','50','60','70','80','90')
	axis([-8 8 -1 3.5])
	grid on
	xlabel('\tau_j')
	ylabel('P/P_{RMS}')
	title('Average Event Shape')
	saveFigure_v2(gcf,[pth 'NoiseEvents-1.5prms\Shape\' F(1:end-4) '_Shape'],300)
	close
	
	
	dtHist = zeros(1,10);
	for n = 1:10
		q = hist(out(n).dt,out(n).dtUnique);
		qi = round(out(n).dtUnique*200000);
		dtHist(qi,n) = q'/sum(q);
	end
	save([pth 'NoiseEvents-1.5prms\' F(1:end-4) '.mat'],'signals','out','D','U','dtHist','thresh')

	figure; hold on;
	CM = colormap('jet');
	for n = 1:10
		plot((1:size(dtHist,1))/200000*U/D,dtHist(:,n),'-o','Color',CM(n*6,:))
	end
	set(gca,'YScale','log')
	legend('90','80','70','60','50','45','40','35','30','25')
	axis([0 4 1e-5 1])
	grid on
	xlabel('\tau_j')
	ylabel('(# Events of Width \tau_j)/(Total Events)')
	title('Histogram of Event Widths')
	saveFigure_v2(gcf,[pth 'NoiseEvents-1.5prms\dtHist\' F(1:end-4) '_dtHist'],300)
	close

	w = signals.w;
	w = reshape(w,8192,100,10);
	w = permute(w,[1 3 2]);
	w = w(:);

	fid = fopen([pth 'NoiseEvents-1.5prms\' F(1:end-4) '.WLT'],'w');
	fwrite(fid,w,'float32');
	fclose(fid);
end
matlabpool close;

save([pth 'NoiseEvents-1.5prms\Nevents.mat'],'Nevents','thresh');
disp('.Done')