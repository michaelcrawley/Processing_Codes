function [f,SPL] = calcSPL(p,Fs,BS,WNDO)
% Calculate the SPL spectrum of p
%
%INPUTS:
% p - a real 2D matrix containing the data arranged with the timeseries
% along the first dimension and multiple channels on the second.
% Fs - Sampling Frequency (Hz)
% BS - Block size of data. Note that the length of the timeseries must be
% an integer multiple of 'BS'.
% WNDO - A appropriate string for the desired windowing of the data. This 
% input is optional - if no input is given, no window is used.

S = size(p);
if round(S(1)/BS)~=S(1)/BS
    error('Data is not an integer number of blocks')
end
p = p - repmat(mean(p,1),S(1),1);   %subtracts mean from each channel
pM = reshape(p,BS,S(1)/BS*S(2));

    %creates window as sparse diagonal matrix
if ~exist('WNDO','var')
    WNDO = 'rectwin';
end
WD = window(str2func(WNDO),BS);
WDM = spdiags(shiftdim(WD),0,BS,BS);
WDW = mean(WD.^2);  %Weight of windowing - must be scaled out of spectrum

nfrq = ceil((BS+1)/2);  %# of distinct frequencies returned for SPL
f = (0:nfrq-1)'*(Fs/BS);   %Frequency axis

rFFT = fft(WDM*pM,BS,1);  %Calculates FFT with the specified window
rFFT = rFFT(1:nfrq,:); %Discards symmetric portion of spectrum
mFFT = abs(rFFT).^2/(BS*Fs*WDW);  %Scales spectrum according to sampling and weighting characteristics 
if rem(BS,2) %Conserves energy which would otherwise be lost by discardiing symmetric portion
    mFFT(2:end,:) = mFFT(2:end,:)*2;
else
    mFFT(2:end-1,:) = mFFT(2:end-1,:)*2;
end
SPL = squeeze(mean(reshape(mFFT,nfrq,S(1)/BS,S(2)),2)); %Calculates squared average of spectrum
