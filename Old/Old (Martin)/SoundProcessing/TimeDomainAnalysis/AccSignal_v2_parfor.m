function [signals,out] = AccSignal_v2_parfor(d,thresh)


% d = load('M13_Baseline_T14.3.NOS');
% CH = 9;
% thresh = cutoff amplitude

sO = [];
TO = [];
wO = [];
signals.ss = [];
signals.T = [];
signals.w = [];

out(size(d,2)).prms = [];
parfor CH = 1:size(d,2)
	dd = d(:,CH);
	
		%This moving average filter extracts ultra-low frequency variation
		%in the signal.
	NM = 500;
	qq = moving_average([ones(NM,1)*mean(dd(1:round(NM)/4)); dd; ones(NM,1)*mean(dd(end-round(NM/4):end))],NM);
	qq = qq(NM+1:length(dd)+NM);
	dd = dd-qq; %clear qq NM;
	
	% [f,SPL,ph] = calcSPL_v2(dd,200000,8192,'hamming');

	% Suppress low-amplitude fluctuations
	S = 10;	%Suppression strength
	prms = std(dd); out(CH).prms = prms;
	I = find(abs(dd) < thresh(CH));
	ss = dd;
	ss(I) = dd(I).*exp(-(thresh(CH)-abs(dd(I)))/(thresh(CH))*S); sO(:,CH) = ss;
	%clear I;

	% Construct Trinary Waveform
	T = zeros(size(dd));
	T(ss > thresh(CH)) = 1;
	T(ss < -thresh(CH)) = -1; TO(:,CH) = T;

	% Find start and end indices for all events
	B = find(diff(abs(T)) > 0);
	E = find(diff(abs(T)) < 0);
	while E(1) < B(1)
		E = E(2:end);
	end
	B = B(1:length(E));

	t = (0:1/200000:(length(dd)-1)/200000)';
	pk = zeros(size(B));    %peak amplitude
	ppms = pk;	%mean square pressure of peak
	for n = length(B):-1:1
		[q,I] = max(abs(dd(B(n):E(n)))); I = I+B(n)-1;
		pk(n) = q*sign(dd(I));
		
		NP = 101;
		q = double(uint64(I-NP))+1;
		tmp = find(dd(q:I)*sign(dd(I)) < pk(n)*sign(dd(I))/2,1,'last');
		while isempty(tmp) && q > 1
			NP = NP+100;
			q = double(uint64(I-NP))+1;
			tmp = find(dd(q:I)*sign(dd(I)) < pk(n)*sign(dd(I))/2,1,'last');
		end
		if isempty(tmp)
			B(n) = 1;
		else
			B(n) = tmp +q-1;
		end
		
		q = I+NP-1;
		if q > length(dd)
			q = length(dd);
		end
		tmp = find(dd(I:q)*sign(dd(I)) < pk(n)*sign(dd(I))/2,1,'first');
		while isempty(tmp) && q < length(dd)
			q = q+100;
			if q > length(dd)
				q = length(dd);
			end
			tmp = find(dd(I:q)*sign(dd(I)) < pk(n)*sign(dd(I))/2,1,'first');
		end
		if isempty(tmp)
			E(n) = length(dd);
		else
			E(n) = tmp +I-1;
		end
		
		ppms(n) = mean(dd(B(n):E(n)).^2)/prms^2;	%Find the average energy contained in a noise event normalized by the average signal energy
		if n < length(B)
			if (E(n) >= B(n+1)) && sign(pk(n))==sign(pk(n+1))
				E(n) = E(n+1);
				[~,I] = max(abs(pk(n:n+1))); pk(n) = pk(n+I-1);

				E = [E(1:n); E(n+2:end)];
				B = [B(1:n); B(n+2:end)];
				pk = [pk(1:n); pk(n+2:end)];
				
				ppms(n) = mean(dd(B(n):E(n)).^2)/prms^2;
				ppms = [ppms(1:n); ppms(n+2:end)];
			end
		end
	end
			
	
	dt = (t(E) - t(B));      %peak width
	to = t(round((B+E)/2));         %peak location
	
	out(CH).dt = dt; out(CH).t0 = to; out(CH).peak = pk; out(CH).ppms = ppms;
	
	[q,~,qb] = unique(E-B); ppms_dtUnique = zeros(size(q)); pss_dtUnique = ppms_dtUnique;
	for n = 1:length(q)	%For each unique width
		ppms_dtUnique(n) = mean(ppms(qb==n));	%Find the average energy of noise events of that width
		pss_dtUnique(n) = sum(ppms(qb==n))*q(n)*prms^2/sum(dd.^2);	%Find the total energy of noise events of that width normalized by the total signal energy
	end
	out(CH).dtUnique = q/200000; out(CH).ppms_dtUnique = ppms_dtUnique; out(CH).pss_dtUnique = pss_dtUnique;

	pw = zeros(size(dd));
	avg = zeros(1000,length(B));
	m = 0;
	for n = 1:length(B)
			%Construct wavelet signal
		w = E(n)-B(n);
		b = B(n)-5*w; b(b < 1) = 1;
		e = E(n)+5*w; e(e > length(pw)) = length(pw);
		pw(b:e) = pw(b:e) +pk(n)*(1 -(t(b:e)-to(n)).^2/(dt(n)*0.8925)^2).*exp(-(t(b:e)-to(n)).^2/(dt(n)*0.8925)^2);
	
	
			%Extract signal for determination of average large event shape
		m = m+1;
		hw = round((E(n)-B(n))*10/2);
		qb = round((B(n)+E(n))/2)-hw;
		qe = round((B(n)+E(n))/2)+hw;
		if qb < 1
			m = m-1;
		elseif qe > length(dd)
			m = m-1;
		else
			if 500-hw < 1
				tmpe = hw-500 +1;
			else
				tmpe = 0;
			end
			avg(500-hw+tmpe:500+hw-tmpe+1,m) = dd(qb+tmpe:qe-tmpe+1)*sign(pk(n));
		end
	end
	wO(:,CH) = pw;
	avg = avg(:,1:m); out(CH).avg = mean(avg,2); out(CH).avgstd = std(avg,[],2);
end

signals.ss = sO; signals.T = TO; signals.w = wO;