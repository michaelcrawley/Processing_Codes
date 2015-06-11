function [f,SPL,ph] = calcSPL_v2(p,FS,BS,WNDO)
% Calculate the SPL spectrum of p.
%
%v2 - 12-10-2009 - improved to include phase and DC information in output.
%This version also includes an inverse transform function.
%
%INPUTS:
% p - a real 2D matrix containing the data arranged with the timeseries
% along the first dimension and multiple channels on the second.
% FS - Hz - Sampling Frequency. If integer indices are desired (i.e.
% integer multiples of 2*pi), set FS = BS.
% BS - Block size of data. Note that the length of the timeseries must be
% an integer multiple of 'BS'. If length(p(:,1))/BS > 1, the function will
% use the multiple blocks to calculate a mean square average spectrum. It
% is strongly recommended that BS be of the form 2^N in order to achieve
% maximum efficiency in "fft.m".
% WNDO - An appropriate string for the desired windowing of the data. This 
% input is optional - if no input is given, no window is used. See the help
% on the function "window" for the list of available windows.
%
%OUTPUTS:
% f - the frequency spectrum for the data
% SPL - amplitude^2 - contains the energy of each frequency in units of
% whatever units p has squared. (e.g. p in units of Pa => SPL in units of
% Pa^2).
% ph - rad - the phase angle for each frequency

S = size(p);
if round(S(1)/BS)~=S(1)/BS
    error('Data is not an integer number of blocks')
end
MN = mean(p,1);
p = p - repmat(MN,S(1),1);   %subtracts mean from each channel - necessary to apply windowing
pM = reshape(p,BS,S(1)/BS*S(2));

    %creates window as sparse diagonal matrix
if ~exist('WNDO','var')
    WNDO = 'rectwin';
end
WD = window(str2func(WNDO),BS);
WDM = spdiags(shiftdim(WD),0,BS,BS);
WDW = mean(WD.^2);  %Weight of windowing - must be scaled out of spectrum

nfrq = ceil((BS+1)/2);  %# of distinct frequencies returned for SPL
f = (0:nfrq-1)'*(FS/BS);   %Frequency axis

rFFT = fft(WDM*pM,BS,1);  %Calculates FFT with the specified window
rFFT = rFFT(1:nfrq,:); %Discards symmetric portion of spectrum
rFFT = reshape(rFFT,nfrq,S(1)/BS,S(2));
mFFT = abs(rFFT).^2/(BS*FS*WDW);  %Scales spectrum according to sampling and weighting characteristics 
if rem(BS,2) %Conserves energy which would otherwise be lost by discardiing symmetric portion
    mFFT(2:end,:,:) = mFFT(2:end,:,:)*2;
else
    mFFT(2:end-1,:,:) = mFFT(2:end-1,:,:)*2;
end
SPL = squeeze(mean(mFFT,2)); %Calculates squared average of spectrum

    %Calculates the phase angle of each frequency. Since phase is not
    %repeatable across non-contiguously sampled blocks, only the first
    %block is used.
tmp = squeeze(rFFT(:,1,:));
ph = angle(tmp); clear tmp;

    %Adds mean information back into spectrum
SPL(1,:) = MN.^2;
ph(1,:) = angle(MN);
