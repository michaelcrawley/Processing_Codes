function [fc,Stc,yc] = thirdOctaveBands(x,St,y)
%Resamples data into 1/3 octave bands.
%
%INPUTS
% x - frequency axis (Hz)
% St - Strouhal axis
% y - SPL in (dB)
%
%OUTPUTS
% fc - 1/3 octave frequencies (Hz)
% Stc - 1/3 octave Strouhal numbers
% yc - 1/3 octave data (dB)

y = 10.^(y/10); %convert to Pa^2

[fc, fc_l, fc_u] = nth_freq_band(3, x(1), x(end));  %Find frequency bands

yc = zeros(size(fc));
for n = 1:length(fc)    %For each frequency band
    IL = min(find(x >= fc_l(n)));   %Locate beginning of band
    IU = max(find(x < fc_u(n)));    %locate end of band
    
    if IL==IU   %If band is only one point wide, use single value.
        yc(n) = y(IL);
    else    %Else the value of the band is the average in the region.
        yc(n) = trapz(x(IL:IU),y(IL:IU))/(x(IU)-x(IL));
    end
end

yc = 10*log10(yc);  %Convert back to dB

Stc = fc*St(1)/x(1);    %Convert Strouhal axis

