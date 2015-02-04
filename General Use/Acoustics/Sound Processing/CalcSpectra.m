function [data] = CalcSpectra(rRaw,pp,data)
%Calculates the spectra from the raw data, and performs any dB corrections
%specified in the processing parameters.

%Last updated by Michael Crawley on 2011-07-04

    [~,SPL] = calcSPL_v2(rRaw,pp.sampleRate,pp.BS,pp.WNDO); 
    
    data.dB = pp.ref+ 10*log10(SPL*pp.ChRef);    %Converts into dB using calibration and reference values
    data.dB = data.dB(pp.DPoints:end-pp.DPoints,:);   %removes junk points at beginning and end
    
    if ~strcmpi('None',pp.atmAbs)
        %Corrects for absorption by atmosphere - Spectra are now representative 
        %of sound propagated to microphone distance with no absorption.
        [alpha,alpha_s] = atmAbsorption_v2(pp.fx,pp.AT,pp.AP,pp.RH,pp.atmAbs); 
        [x,y] = meshgrid(pp.ChD,alpha); data.dB = data.dB + x.*y;    %Removes absorption at current conditions
        clear alpha x y;
    end
    
    if strcmpi('Radial',pp.NDs)    %Correct to normalized distance using simple radial propagation
        if sum(abs(pp.MO))~=0
            data.dBCorrected = micOrientationCorrect(pp.fx,pp.MO,data.dB); %Corrects spectra for microphone orientation
        else
            data.dBCorrected = data.dB;
        end
        [data.dBCorrected, al] = NormDistance(data.dBCorrected,pp.ChD/pp.D,pp.NormD,pp.dBShift);
        if strcmpi('Standard',pp.atmAbs)   %Adds absorption at standard day conditions if selected
            [x,y] = meshgrid(pp.ChD.*al',alpha_s); data.dBCorrected = data.dBCorrected - x.*y;
        end
    elseif strcmpi('NS',pp.NDs)        %Correct to normalized distance using noise source location

        [data.dBCorrected, al, ph] = NormDistance(data.dB,data.Std,pp.ChD/pp.D,pp.ChPol*pi/180,pp.NormD,pp.dBShift);
        if sum(abs(pp.MO))~=0
            MOph = abs(ph*180/pi - repmat(pp.MO+pp.ChPol,length(data.Std),1));
            data.dBCorrected = micOrientationCorrect(repmat(pp.fx,1,length(pp.MO)),MOph,data.dBCorrected); %Corrects spectra for microphone orientation        
        end
        if strcmpi('Standard',pp.atmAbs)
            [x,y] = meshgrid(pp.ChD,alpha_s); data.dBCorrected = data.dBCorrected - x.*y.*al;
        end
    else
        if sum(abs(pp.MO))~=0
            data.dBCorrected = micOrientationCorrect(pp.fx,pp.MO,data.dB);
        else
            data.dBCorrected = data.dB;
        end
        if strcmpi('Standard',pp.atmAbs)   %Adds absorption at standard day conditions if selected
            [x,y] = meshgrid(pp.ChD,alpha_s); data.dBCorrected = data.dBCorrected - x.*y;
        end
    end
    
    data.dBCorrected_DT = data.dBCorrected; 
    %Smoothes and removes forcing harmonics from spectrum. If FH==0 (forcing frequency) program doesn't look for harmonics
    if pp.BT ~= 0 %if the tones are broadband (pp.BT is true), use alternate tone subtract method
        if pp.DT ~= 0
            for m = 1:pp.Nch
                data.dBCorrected_DT(:,m) = ToneSubtract_v5(pp.fx,data.dBCorrected(:,m),'Moving',data.FH,10);
            end
        else
            for m = 1:pp.Nch
                data.dBCorrected_DT(:,m) = ToneSubtract_v5(pp.fx,data.dBCorrected(:,m),'Moving',0,10);
            end
        end        
    else %if the tones are sharp (pp.BT is false), use standard tone subtract method
        for m = 1:pp.Nch   
            data.dBCorrected_DT(:,m) = ToneSubtract_v5(pp.fx,data.dBCorrected(:,m),'Poly',30);
        end
    end
end