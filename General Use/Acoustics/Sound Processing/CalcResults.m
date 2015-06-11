function [results] = CalcResults(data,pp) 
%Calculates OASPL, AAE, NJP for given data set.  If delta OASPL, delta AAE etc are
%desired, baseline results need to be included in the processing parameters
%structure in the same format as the results structure, named as 'bdata'.
%This can be accomplished by running the ProcessBaselineResults.m function.

%Last updated by Michael Crawley on 2011-07-04

    %Calculates OASPL for all channels    
    Bg = find(data.Std > pp.iLimits(1),1,'first'); %Excludes spectral components below 0.04 Strouhal number
    Ed = find(data.Std > pp.iLimits(2),1,'first');    %Excludes spectral components above 4 Strouhal number
    if isempty(Bg)
        Bg = 1;
    end
    if isempty(Ed)
        Ed = length(data.Std);
    end
    results.OASPL = 10*log10(trapz(pp.fx(Bg:Ed),10.^(data.dBCorrected(Bg:Ed,:)/10),1)); 
    
        %Calculates the average acoustic energy
    delta_Pol = abs(diff(pp.ChPol)); IW = ([delta_Pol/2 0] + [0 delta_Pol/2])/(max(pp.ChPol)-min(pp.ChPol));
    results.AAE = 10*log10(sum(IW.*10.^(results.OASPL/10)));
   
        %Calculates the normalized jet total power - assumes
        %azimuthally symmetric and fully expanded. The solid angle of the 
        %integration is divided out giving (W/sr). This power per solid
        %angle is then normalized by the jet mechanical power per solid
        %angle of a complete sphere (i.e. 4pi).
    results.NJP = 10*log10(32*pp.NormD^2*pp.pref^2*sqrt((pp.AT+273.15)/data.Temp)*...
        trapz(pp.ChPol*pi/180,10.^(results.OASPL/10).*sin(pp.ChPol*pi/180))/...
        ((1.4*pp.AP*1000)^2*pp.M^3*(cos(pp.ChPol(1)*pi/180)-cos(pp.ChPol(end)*pi/180))));
        

    iPeakSPL = zeros(1,pp.Nch);
    results.PeakSPL = iPeakSPL;
    results.PeakSTD = iPeakSPL;
    results.PE = iPeakSPL;
    results.Pval = zeros(2,pp.Nch);
    I = find(data.Std >= 0.5,1,'first');   %Ignores all points to left of Std = 0.5
    for nn = 1:pp.Nch
       %Calculates the peak SPL magnitude and Strouhal number
       [~,iPeakSPL(nn)] = max(data.dBCorrected_DT(1:Ed,nn));
       results.PeakSTD(nn) = data.Std(iPeakSPL(nn));
       iPeakSPL(iPeakSPL < 6) = 6;
       results.PeakSPL(nn) = mean(data.dBCorrected_DT(iPeakSPL(nn)-5:iPeakSPL(nn)+5,nn));
       
       %Calculates the high frequency slope of the spectra including uncertainty as seen on semilogx scale.
       results.Pval(:,nn) = mmpolyfit(log10(data.Std(I:Ed)),data.dBCorrected_DT(I:Ed,nn),1,'Weight',1./data.Std(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
       %Calculates uncertainty of fit
       Serr = sum((results.Pval(1,nn)*log10(data.Std(I:Ed))+results.Pval(2,nn)-data.dBCorrected_DT(I:Ed,nn)).^2);
       merr = length(data.Std(I:Ed));
       Derr = merr*sum(log10(data.Std(I:Ed)).^2) - sum(log10(data.Std(I:Ed)))^2;
       results.PE(nn) = sqrt(Serr*merr/(merr-2)/Derr);        
    end
    
    %If forcing is present, calculates/prints several relevant values
    if data.SHarmonics
        %Calculates OASPL for all channels
        if pp.DT~=0
            results.OASPL_DT = 10*log10(trapz(pp.fx(Bg:Ed),10.^(data.dBCorrected_DT(Bg:Ed,:)/10),1)); %clear dBCorrected_DT; [Modification by M Crawley]
                %Calculates the average acoustic energy
            results.AAE_DT = 10*log10(sum(IW.*10.^(results.OASPL_DT/10)));
                %Calculates and prints the normalized jet total power
            results.NJP_DT = 10*log10(32*pp.NormD^2*pp.pref^2*sqrt((pp.AT+273.15)/data.Temp)*...
                trapz(pp.ChPol*pi/180,10.^(results.OASPL_DT/10).*sin(pp.ChPol*pi/180))/...
                ((1.4*pp.AP*1000)^2*pp.M^3*(cos(pp.ChPol(1)*pi/180)-cos(pp.ChPol(end)*pi/180))));
        end
        
        if isfield(pp,'bdata')                
            if length(pp.bdata.BT) > 1
                for nn = 1:pp.Nch      %Interpolate baseline results to most exactly match baseline characteristics at "To"
                    pp.bdata.OASPLB(1,nn) = interp1(pp.bdata.BT,pp.bdata.OASPLB(:,nn),data.To,'linear','extrap');
                    pp.bdata.OASPLB_DT(1,nn) = interp1(pp.bdata.BT,pp.bdata.OASPLB_DT(:,nn),data.To,'linear','extrap');
                end
                pp.bdata.OASPLB = pp.bdata.OASPLB(1,:);
                pp.bdata.OASPLB_DT = pp.bdata.OASPLB_DT(1,:);
                pp.bdata.AAEB = interp1(pp.bdata.BT,pp.bdata.AAEB,data.To,'linear','extrap');
                pp.bdata.AAEB_DT = interp1(pp.bdata.BT,pp.bdata.AAEB_DT,data.To,'linear','extrap');
                pp.bdata.NJPB = interp1(pp.bdata.BT,pp.bdata.NJPB,data.To,'linear','extrap');
                pp.bdata.NJPB_DT = interp1(pp.bdata.BT,pp.bdata.NJPB_DT,data.To,'linear','extrap');
            end
            
                            %Calculates delta OASPL
            results.dOASPL = results.OASPL-pp.bdata.OASPLB;
                %Calculates change in average acoustic energy
            results.dAAE = results.AAE-pp.bdata.AAEB;
                %Calculates change in normalized jet power
            results.dNJP = results.NJP-pp.bdata.NJPB;

            if pp.DT~=0
                    %Calculates delta OASPL
                results.dOASPL_DT = results.OASPL_DT - pp.bdata.OASPLB_DT;
                    %Calculates change in average acoustic energy
                results.dAAE_DT = results.AAE_DT-pp.bdata.AAEB_DT;
                    %Calculates change in normalized jet power
                results.dNJP_DT = results.NJP_DT-pp.bdata.NJPB_DT;
            end 
        end    
    else
            %Calculates OASPL for all channels    
        results.OASPL_DT = 10*log10(trapz(pp.fx(Bg:Ed),10.^(data.dBCorrected_DT(Bg:Ed,:)/10),1)); 
            %Calculates the average acoustic energy
        results.AAE_DT = 10*log10(sum(IW.*10.^(results.OASPL_DT/10)));
            %Calculates the normalized jet total power
        results.NJP_DT = 10*log10(32*pp.NormD^2*pp.pref^2*sqrt((pp.AT+273.15)/data.Temp)*...
            trapz(pp.ChPol*pi/180,10.^(results.OASPL_DT/10).*sin(pp.ChPol*pi/180))/...
            ((1.4*pp.AP*1000)^2*pp.M^3*(cos(pp.ChPol(1)*pi/180)-cos(pp.ChPol(end)*pi/180))));
    end
end