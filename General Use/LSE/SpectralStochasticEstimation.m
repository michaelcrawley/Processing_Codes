function[cond_coeff,fieldr,nr,nc,nrs,ncs] = SpectralStochasticEstimation(signal,fss,field,fsf,nblock,zero_padding)

% 'signal' and 'field' should be reshaped as a 2D field with the second
% dimension (n column) as the time.

if nargin < 5
    nblock = 10;
    zero_padding = 0;
elseif nargin < 6
    zero_padding = 0;
end
[nr,nc] = size(field);
[nrs,ncs] = size(signal);
signal = signal - repmat(mean(signal,2),1,ncs);
% Down-sampling signal
signaltmp = signal(:,1:fss/fsf:(fss/fsf)*nc); clear signal
signal = signaltmp; clear signaltmp
ns = floor(size(signal,2)/nblock); % Size of each block

% Taking out the mean value of the field
% field = field - repmat(mean(field,2),1,nc);

% Divinding 'signal' and 'field' in blocks
signalb = zeros(nrs,ns,nblock);
fieldb = zeros(nr,ns,nblock);
for it_b = 1:nblock
    signalb(:,:,it_b) = signal(:,1+(it_b-1)*ns:ns*it_b);
    fieldb(:,:,it_b) = field(:,1+(it_b-1)*ns:ns*it_b);
end, clear it_b

% FFT of 'signal' and 'field', with or without 'zero padding'
if zero_padding == 0
    fsignal = fft(signalb,[],2);
    ffield = fft(fieldb,[],2);
else
    fsignal = fft(signalb,2*ns,2);
    ffield = fft(fieldb,2*ns,2);
end

unconditional = complex(zeros(nrs,nrs,size(fsignal,2)));
for it_frq = 1:size(fsignal,2)
    for it_b = 1:nblock
        unconditional(:,:,it_frq) = unconditional(:,:,it_frq) + fsignal(:,it_frq,it_b)*(fsignal(:,it_frq,it_b)'); % Matrix W
    end, clear it_b
    unconditional(:,:,it_frq) = unconditional(:,:,it_frq)/nblock;
end, clear it_frq

conditional = complex(zeros(nr,nrs,size(fsignal,2)));
for it_frq = 1:size(fsignal,2)
    for it_b = 1:nblock
        conditional(:,:,it_frq) = conditional(:,:,it_frq) + ffield(:,it_frq,it_b)*(fsignal(:,it_frq,it_b)'); % Matrix V
    end, clear it_b
    conditional(:,:,it_frq) = conditional(:,:,it_frq)/nblock;
end, clear it_frq

cond_coeff = complex(zeros(nr,size(unconditional,2),size(fsignal,2)));
for it_frq = 1:size(fsignal,2)
    cond_coeff(:,:,it_frq) = (unconditional(:,:,it_frq)\(conditional(:,:,it_frq).')).'; % Matrix A
end, clear it_frq

ffieldr = complex(zeros(nr,size(fsignal,2),nblock));
for it_frq = 1:size(fsignal,2)
    for it_b = 1:nblock
        ffieldr(:,it_frq,it_b) = cond_coeff(:,:,it_frq)*(fsignal(:,it_frq,it_b));
    end, clear it_b
end, clear it_frq

fieldr = zeros(nr,ns*nblock);
for it_b = 1:nblock
    fieldr(:,1+(it_b-1)*ns:it_b*ns) = ifft(ffieldr(:,:,it_b),[],2);
end, clear it_b