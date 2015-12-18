function [t,p] = icalcSPL_v2(SPL,ph,FS,BS,WNDO)
% Calculate the data p from SPL spectrum. Note that this inverse transform
% is not capable of "undoing" the RMS averaging over multiple blocks. In
% the case of multiple blocks used in the forward transform, the output 
% will be 1 block which is, by definition, the data corresponding to the
% RMS averaged spectrum.
%
%v2 - 12-10-2009 - First release of inverse transform function.
% -NOTE - At this point, the inversion function does not properly handle
% windows. It is recommended that this function only be used with no
% window.
%
%INPUTS:
% SPL - amplitude^2 - contains the energy of each frequency in units of
% whatever units p has squared. (e.g. p in units of Pa => SPL in units of
% Pa^2).
% ph - rad - the phase angle for each frequency
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
%----NOTE: FS, BS, and, WNDO must be the same as those used in the forward
%transform in order to obtain the correct inversion.
%
%OUTPUTS:
% t - the frequency spectrum for the data
% p - a real 2D matrix containing the data arranged with the timeseries
% along the first dimension and multiple channels on the second.

%% THIS CODE DOESN"T WORK %%%%%%%%%%%%%%%%%%%%
    %creates window as sparse diagonal matrix
% if ~exist('WNDO','var')
%     WNDO = 'rectwin';
% end
% WD = 1./window(str2func(WNDO),BS);  %Creates inverted window to remove windowing effect from data
% WDM = spdiags(shiftdim(WD),0,BS,BS);
% WDW = 1/mean(WD.^2);  %Weight of windowing - must be scaled out of spectrum
WDW = 1;

%%
t = (0:1/FS:(BS-1)/FS)';    %Calculates real-space coordinates of signal

MN = real(sqrt(SPL(1,:)).*exp(i*ph(1,:)));    %Extracts mean information from spectrum
SPL(1,:) = 0;   %Removes mean from spectrum - necessary to apply windowing

mFFT = SPL*(BS*FS*WDW); %Removes spectrum scaling according to sampling and weighting characteristics
if rem(BS,2) 
    mFFT(2:end,:) = mFFT(2:end,:)/2;    %Conserves energy which would otherwise be lost by discardiing symmetric portion
    rFFT = [sqrt(mFFT).*exp(i*ph); flipud(sqrt(mFFT(2:end,:)).*exp(-i*ph(2:end,:)))];   %Unfolds the spectrum into the form required by "ifft.m"
else
    mFFT(2:end-1,:) = mFFT(2:end-1,:)/2;
    rFFT = [sqrt(mFFT).*exp(i*ph); flipud(sqrt(mFFT(2:end-1,:)).*exp(-i*ph(2:end-1,:)))];
end

iFT = ifft(rFFT,BS,1);  %Computes the iFFT
% p = WDM*real(iFT) + repmat(MN,BS,1);  %Removes the window and adds the DC component
p = real(iFT) + repmat(MN,BS,1);  %adds the DC component
