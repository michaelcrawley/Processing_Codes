
tic;
pool = parpool(12);
%Preprocessing for wavelet transform
smoothed2 = p_s;
% signal = p_s(1:128,9:520,:);
qmf = MakeONFilter('Battle',5);
[M,N,L] = size(smoothed2);

% %Filter in radial direction
% signal = permute(signal,[2 3 1]);
% signal = reshape(signal,N*L,M);
% parfor k = 1:(N*L)
%     tmp = signal(k,:);
%     filt = CohWave(tmp,4,qmf);
%     signal(k,:) = filt;
% end
% signal = reshape(signal,[N,L,M]);
% signal = ipermute(signal,[2 3 1]);
% 
% %Filter in axial direction
% signal = permute(signal,[1 3 2]);
% signal = reshape(signal,M*L,N);
% parfor k = 1:(M*L)
%     tmp = signal(k,:);
%     filt = CohWave(tmp,4,qmf);
%     signal(k,:) = filt;
% end
% signal = reshape(signal,[M,L,N]);
% signal = ipermute(signal,[1 3 2]);

% Smooth in spatial dimensions
disp('Smoothing in Space');
parfor k = 1:L
    smoothed2(:,:,k) = moving_average2(smoothed2(:,:,k),1,2);
end

% 
% %Filter in temporal direction
% disp('Smoothing in Time');
% smoothed1 = reshape(smoothed1,[],L);
% smoothed2 = reshape(smoothed2,[],L);
% parfor k = 1:(M*N)
%     tmp = smoothed1(k,:);
%     filt = CohWave(tmp,4,qmf);
%     smoothed1(k,:) = filt;
%     
%     tmp = smoothed2(k,:);
%     filt = ThreshWave(tmp,'S',1,std(tmp),multi,2,qmf);
%     smoothed2(k,:) = filt;
% end
% smoothed1 = reshape(smoothed1,[M,N,L]);
% smoothed2 = reshape(smoothed2,[M,N,L]);

% % Smooth in spatial dimensions
% disp('Smoothing in Space');
% multi = sqrt(log(2^max(nextpow2([M N]))^2));
% parfor k = 1:L
%     tmp = zeros(2^max(nextpow2([M N])));
%     tmp(1:M,1:N) = smoothed(:,:,k);
%     tmp = ThreshWave2(tmp,'S',1,std(tmp(:)),multi,3,qmf);    
%     smoothed(:,:,k) = tmp(1:M,1:N);
% end

%Filter in temporal direction
disp('Smoothing in Time');
multi = 1.5;%sqrt(log(L)-1);
smoothed2 = reshape(smoothed2,[],L);
parfor k = 1:(M*N)    
    tmp = smoothed2(k,:);
    filt = ThreshWave(tmp,'S',1,std(tmp),multi,2,qmf);
    smoothed2(k,:) = filt;
end
smoothed2 = reshape(smoothed2,[M,N,L]);


delete(pool);
toc