function [SPL, f] = calSPL(p, nfft, fs, WindStrFnc, noverlap, flags)
% Calculate the SPL spectrum of p

%   Assumptions:
%   1.  'p' is provided in Pascals.
%   2.  'p' is a 2D matrix with timeseries along 1st dim and channels along
%       2nd dim.
%   3.  'p' can be complex (as when it is an azimuthally Fourier decomposed
%       signal), in which case p(-f) ~= p*(f) in general. But typical
%       symmetries are such that, we still expect SPL(-f) = SPL(f). So, in 
%       all cases, we will return only the SPL over the +f axis, and double
%       the energy in it to account for the neglected -f axis.
%   4.  For overlapping, it is of course assumed that the timeseries is
%       contiguous.

%   Input Parameters:
%   p:          time-series 
%   nfft:       # frequency points used to calculate the discrete Fourier
%               transforms
%   fs:         sampling frequency, in Hz
%   WindStrFnc: String indicating window function to use. If it is not
%               provided, or it is an empty string, then it defaults to the
%               rectangular window. Otherwise, the string is converted to
%               the function handle and used for generating the window.
%   noverlap:   # samples each segment of p overlaps. Defaults to 0.
%   flags:      flags(1):   1:  Normalize by p_ref (default)
%                           0:  Don't normalize
%               flags(2):   1:  Give result in dB
%                           0:  Don't convert to dB

%   Output Parameters:
%   SPL:        Sound Pressure Level, in dB.
%               It is a 2D matrix with 'nfrq' (see below) frequencies along
%               1st dim and channels along 2nd dim.
%   f:          nfrq long vector of frequency axis

%   NOTE:
%   The same variable is returning different quantities depending on the 
%   setting of "flags"; hence its name is ambiguous. Use the "flags" input
%   to determine the correct meaning. Eg. when flags(1) = 0, the output means
%   PSD.

%   If WindStrFnc is not provided, or empty, then use the rectangular window
if nargin < 4 || isempty(WindStrFnc)
    WindStrFnc = 'rectwin';
end
%   Window function
dat_wind = window(str2func(WindStrFnc),nfft);
dat_wind_mtx = spdiags(shiftdim(dat_wind),0,nfft,nfft);
wind_weight = mean(dat_wind.^2);

%   If noverlap is not provided, or empty, then defaults to 0
if nargin < 5 || isempty(noverlap)
    noverlap = 0;
end

if nargin < 6 || isempty(flags)
    flags = [1,1];
elseif length(flags) == 1
    flags = [flags 1];
end

if flags(1)
    %   Reference pressure (20 micro-Pascal)
    pref = 2e-5;
    %   Normalized pressure timetrace
    p = p/pref;
end

%   # distinct frequencies returned for SPL
nfrq = ceil((nfft+1)/2);

%   Frequency axis
f = (0:nfrq-1)*(fs/nfft);

%   # channels
nch = size(p,2);
%   # samples by which one block is shifted from the next
blk_shift = nfft - noverlap;
%   # blocks to perform averaging over
nblks = fix((size(p,1)-noverlap)/blk_shift);

%   Pre-allocate PSD matrix for efficiency. Note that we allocate for nfft
%   frequencies initially. We will do the fold-over subsequently.
PSD = zeros(nfft,nch);
% Perform calculations separately for each channel
for blk = 1:nblks
    blk_curr = (1:nfft)+(blk-1)*blk_shift;
    p_f = fft(dat_wind_mtx*p(blk_curr,:),nfft,1);
    PSD = PSD + (abs(p_f)).^2;
end
PSD = PSD/nblks;
%   Account for the neglected energy. For real signals PSD(-f) = PSD(f).
%   But this is not true for complex signals in general. However, we assume
%   that it is approximately true (which is true in all cases that we
%   typically consider), so that a direct summation is warranted.
if rem(nfft,2)
    PSD(2:nfrq,:) = PSD(2:nfrq,:) + PSD(end:-1:nfrq+1,:);
else
    PSD(2:nfrq-1,:) = PSD(2:nfrq-1,:) + PSD(end:-1:nfrq+1,:);
end

SPL = PSD(1:nfrq,:)/(nfft*fs*wind_weight);
if flags(2)
    %   SPL, in dB, normalized appropriately.
    SPL = 10*log10(SPL);
end
