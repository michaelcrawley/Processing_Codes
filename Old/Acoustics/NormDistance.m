function [Aout, al] = NormDistance(varargin)
%Corrects acoustic results for distance based on noise source location or
%simple radial propagation depending on inputs. The noise source location
%calculation is at the end of this file.
%
%CALLS
% [Aout, al] = NormDistance(A,Std,R,ChPol,NormD,dBShift)
%       - Corrects based on noise source location
% [Aout, al] = NormDistance(A,R,NormD,dBShift)
%       - Corrects using simple radial propagation
%
%INPUTS
% A - Amplitude (db) - Should be m-rows by n-channels
% Std - Strouhal number 
% R - vector containing radial distance to n-microphones
% ChPol - polar angle from jet axis to microphone (rad)
% NormD - Normalization distance 
% dBShift - arbitrary uniform offset shift (dB)
%
%OUTPUTS
% Aout - corrected data (dB)
% al - distance multiplier for propogated distance

switch nargin
    case 4
        A = varargin{1};
        R = varargin{2};
        NormD = varargin{3};
        dBShift = varargin{4};
    case 6
        A = varargin{1};
        Std = varargin{2};
        R = varargin{3};
        ChPol = varargin{4};
        NormD = varargin{5};
        dBShift = varargin{6};
    otherwise
        error('Bad inputs')
end

if nargin==4
    al = NormD./R(:);
    ChDM = spdiags(20*log10(1./al) +dBShift,0,length(R),length(R)); %prepares normalization distance for scaling the data.
    Aout = A+ ones(size(A))*ChDM;     %Normalizes data to a standard distance
else
    x = NoiseSourceDist(Std); %find locations for noise sources by frequency
    [l,Rd] = pol2cart(ChPol(:)',R(:)');
    % l = Rd./tan(ChPol);
    % R = Rd./sin(ChPol);
    [X L] = meshgrid(x,l);
    s = abs(X-L);
    CR = sqrt(repmat(Rd(:).^2,[1 length(x)])+s.^2); %determine corrected r/D for each microphone & frequency pair
    CP = atan(repmat(Rd(:),[1 length(x)])./s); %determine corrected polar angle for each microphone & frequency pair
    tempP = meshgrid(ChPol,1:length(Std))';
    tempR = meshgrid(R,1:length(Std))';
    beta = (pi+tempP-CP).*(X<=L)+(tempP+CP).*(X>L);
    dx = tempR.*cos(beta)+sqrt((tempR.*cos(beta)).^2-tempR.^2+NormD^2); %determine added distance to normalize to, from NSL
    CND = CR+dx;%determines normalization distance for each frequency/microphone cell    
    ph = (acos((X+CND.*cos(CP))/NormD).*(X < L)+acos((X-CND.*cos(CP))/NormD).*(X > L))';%determine angle between normalization (from NSL) and nozzle exit
    al = (CND./CR)';
    At = A + 20*log10(1./al)+dBShift;
    Aout = zeros(size(A));   

    for m = 1:length(Std)    %Determine correct propagation for each frequency since source location is a function of frequency   
        for n = 1:length(ChPol)     %Uses the calculated angle data "ph" to determine a weighted average for accurate propagation
            pA = find(ph(m,:) > ChPol(n),1,'last');
            pB = find(ph(m,:) <= ChPol(n),1,'first');
            if ~isempty(pA) && ~isempty(pB)
                dph = ph(m,pA)-ph(m,pB);
                W = 1 -(ph(m,pA)-ChPol(n))/dph;
                Aout(m,n) = W*At(m,pA) + (1-W)*At(m,pB);
            elseif ~isempty(pA)
                Aout(m,n) = At(m,pA);
            else
                Aout(m,n) = At(m,pB);
            end
        end
    end
end

return



function x = NoiseSourceDist(Std)
    xD = @(Std) 5.1-7.3*log10(Std); %define function to determine source location for strouhal number, from Kuo 2010
    x = xD(Std); %find locations for noise sources by frequency
    x(x<0) = 0;
return