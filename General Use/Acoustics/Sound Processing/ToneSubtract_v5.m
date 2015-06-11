function yo = ToneSubtract_v5(x,y,SW,varargin)
%
%Version 5 adds switch (SW) for choice of polynomial fitting. Polynomial
%fitting works better for spectra with narrowband tones. If tones are more
%broadband, MovingAverage subtraction is better.
%
%Removes harmonics of forcing frequency and smoothes spectrum for ease of
%visualization and calculation of OASPL. This algorithm implements a
%variable size window to better smooth spectra commonly viewed on a semilog
%x-axis given by: 10*log10(x/x(N+1))+N. This is a multiple pass
%implementation using both forward and backward passes to avoid shifting
%peak values. Additionally, harmonics of the forcing frequency (assumed to
%be a square wave with 10% duty cycle) are determined and removed by local
%point negation so that the averaging eliminates the feature instead of
%flattening and broadening it.
%
%Set F=0 to eliminate harmonic removal segment of program.
%
% INPUTS
%x - Horizontal axis (Hz) - should be frequency (not Strouhal number)
%y - SPL - the spectrum - single vector only
%SW - switch for choice of polynomial ('Poly') or moving average ('Moving')
% smoothing technique.
%pO - polynomial order
%F - Forcing frequency (Hz) - should be the frequency set in the
% controller. This algorithm corrects for slight frequency discrepancy of
%1.5% observed in controller operation.
%N - side band size - the number of points on either side of the point
% being averaged at low frequency limit

% N = 10;   %typical value for N
%
%Function calls:
% yo = ToneSubtract_v5(x,y,'Poly',pO);
% yo = ToneSubtract_v5(x,y,'Moving',F,N);

if strcmpi(SW,'Poly')
    SW = 1;
    pO = varargin{1};
else
    SW = 0;
    F = varargin{1};
    N = varargin{2};
end

if SW
    warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
    Cp = polyfit(x,y,pO);
    yo = polyval(Cp,x);
else
        %Harmonic inversion segment
    if F~=0
        fm = (1:1:max(x)/F)*F;
        fm = fm(fm < x(end)*0.99);  %skips harmonic inversion for tones near upper limit of spectrum
        if F >= 5000   %For forcing frequencies above 5 kHz, keeps at most first 7 harmonics and performs harmonic inversion
            fm = fm(fm < 7*F);  %throws out harmonics above 7 times fundamental
            for n = 1:length(fm) 
                fi = min(find(x > fm(n)));
                q = mean(y(fi-50:fi-20));
                y(fi-2:fi+2) = 2*q - y(fi-2:fi+2);
            end
        else
            fm = fm(fm < 15*F);  %Keeps at most 15 harmonics
            for n = 1:length(fm)
                fi = min(find(x > fm(n)));
                if fm(n) < 2000 %performs harmonic inversion using different window sizes based on frequency of harmonic
                    q = (mean(y(fi-10:fi-5))+mean(y(fi+5:fi+10)))/2;
                    y(fi) = 2*q - y(fi);
                elseif fm(n) < 5000
                    q = mean(y(fi-50:fi-20));
                    y(fi-1:fi+1) = 2*q - y(fi-1:fi+1);
                else
                    q = mean(y(fi-50:fi-20));
                    y(fi-3:fi+3) = 2*q - y(fi-3:fi+3);
                end
            end
        end
    end

        %Low pass filtering segment
    NV = round(10*log10(x/x(N+1))+N);   %Generates variable window size 
    yo = y;
    for n = N+1:length(y)-NV(end)   %forward pass
        yo(n) = mean(yo(n-NV(n):n+NV(n)));
    end
    for n = length(y)-NV(end):-1:N+1    %backward pass
        yo(n) = mean(yo(n-NV(n):n+NV(n)));
    end

    if F~=0 %additional two passes on frequencies above 10 kHz if harmonic inversion was done
        S = min(find(x > 10000));
        NV = round(NV.*(log10(x/x(S))+1));  %uses enlarged window size
        S = min(find(NV > 0));  %ensures second set of passes blends smoothly with piece of spectrum which is not additionally smoothed
        for n = S:length(y)-NV(end) %forward pass
            yo(n) = mean(yo(n-NV(n):n+NV(n)));
        end
        for n = length(y)-NV(end):-1:S  %backward pass
            yo(n) = mean(yo(n-NV(n):n+NV(n)));
        end
    end

    for n = N:-1:1  %smoothes low end of spectrum
        yo(n) = mean(yo(1:2*n));
    end

    for n = length(y)-NV(end)+1:length(y)   %smoothes high end of spectrum
        yo(n) = mean(yo(2*n-length(y):end));
    end
end
