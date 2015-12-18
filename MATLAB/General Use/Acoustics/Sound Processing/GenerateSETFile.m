function GenerateSETFile(filename,results,data,pp,set_dir) 
%Generates *.SET files based on values in 'pp', 'data', and
%'results', structures.  *.SET files are saved in the directory specified
%by set_dir.

%Last updated by Michael Crawley on 2011-07-04

    fid = fopen([set_dir '\' filename(1:end-3) 'SET'],'w');
    fprintf(fid,'PROCESS VALUES FOR DATA\n');
    fprintf(fid,['Processing Time Stamp:' pp.delimO datestr(now,0) '\n']);
    fprintf(fid,['Sampling Rate (Hz):' pp.delimO num2str(pp.sampleRate) '\n']);
    fprintf(fid,['Block Size:' pp.delimO num2str(pp.BS) '\n']);
    fprintf(fid,['Reference Amplitude (dB):' pp.delimO num2str(pp.ref) '\n']);
    fprintf(fid,['Window:' pp.delimO pp.WNDO '\n']);
    fprintf(fid,['Normalization Distance (x/D):' pp.delimO num2str(pp.NormD) '\n']);    
    fprintf(fid,['Normalization Method:' pp.delimO pp.NDs '\n']);
    fprintf(fid,['Number of Channels:' pp.delimO num2str(pp.Nch) '\n']);
        %Records location and name of each calibration file used.
    fprintf(fid,['Location of Calibration Files:' pp.delimO pp.CalPathP '\n']);
    fprintf(fid,'Calibration Files:');
    for nn = 1:pp.Nch
        fprintf(fid,[pp.delimO pp.CalFiles{nn}]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mic Distance (m):'); fprintf(fid,'\t%.2f',pp.ChD); fprintf(fid,'\n');
    fprintf(fid,'Mic Polar Angle (deg):'); fprintf(fid,'\t%.2f',pp.ChPol); fprintf(fid,'\n');
    fprintf(fid,'Mic Orientation Angle (deg):'); fprintf(fid,'\t%.2f',pp.MO); fprintf(fid,'\n');
    if pp.DT==0
        fprintf(fid,'Tone Subtraction:\tNO\n');
    else
        fprintf(fid,'Tone Subtraction:\tYES\n');
    end
    fprintf(fid,['Spectrum offset (dB):' pp.delimO num2str(pp.dBShift) '\n']);
    fprintf(fid,['Trim Ends of Spectrum (N):' pp.delimO num2str(pp.DPoints) '\n']);
    if strcmpi('Standard',pp.atmAbs)
        fprintf(fid,'Atmospheric Absorption:\tStandard Day\n');
    elseif strcmpi('Lossless',pp.atmAbs)
        fprintf(fid,'Atmospheric Absorption:\tLossless Propagation\n');
    else
        fprintf(fid,'Atmospheric Absorption:\tNo Correction\n');
    end
    fprintf(fid,['Jet Diameter (m):' pp.delimO num2str(pp.D) '\n']);
    fprintf(fid,['Mach Number:' pp.delimO num2str(pp.M) '\n']);
    fprintf(fid,['Ambient Temperature (K):' pp.delimO num2str(pp.AT+273.15) '\n']);
    fprintf(fid,['Ambient Pressure (kPa):' pp.delimO num2str(pp.AP) '\n']);
    fprintf(fid,['Relative Humidity (%%):' pp.delimO num2str(pp.RH) '\n']);
    fprintf(fid,['Stagnation Temperature (K):' pp.delimO num2str(data.To) '\n']);
    fprintf(fid,['Exit Temperature (K):' pp.delimO num2str(data.Temp) '\n']);
    fprintf(fid,['Total Temperature Ratio:' pp.delimO num2str(data.TTR) '\n']);
    fprintf(fid,['Exit Velocity (m/s):' pp.delimO num2str(data.Ue) '\n']);
    fprintf(fid,['Jet Velocity Ratio (Uj/a_inf):' pp.delimO num2str(data.Mj) '\n']);
    fprintf(fid,['Exit Density (kg/m^3):' pp.delimO num2str(data.rhoe) '\n']);
    fprintf(fid,['Exit Viscosity (Pa s):' pp.delimO num2str(data.mue) '\n']);   %equation is Sutherlands formula - 291.15 is reference temperature
    fprintf(fid,['Jet Reynolds Number:' pp.delimO num2str(data.ReN,'%.0f') '\n']);
    
    fprintf(fid,'OASPL (dB):'); fprintf(fid,'\t%.3f',results.OASPL); fprintf(fid,'\n');
    fprintf(fid,['Average Energy (dB):' pp.delimO num2str(results.AAE,'%.3f') '\n']);
    fprintf(fid,['Normalized Jet Power (dB):' pp.delimO num2str(results.NJP,'%.3f') '\n']);
      
    fprintf(fid,'Peak of SPL Curve\n');    
    fprintf(fid,'  Magnitude (dB):');
    fprintf(fid,[pp.delimO '%.3f'],results.PeakSPL);
    fprintf(fid,'\n');
    
    fprintf(fid,'  Frequency (Std):');
    fprintf(fid,[pp.delimO '%.4f'],results.PeakSTD);
    fprintf(fid,'\n');
    
    fprintf(fid,'Slope of SPL Curve (Strouhal axis):');
    for nn = 1:pp.Nch
        fprintf(fid,[pp.delimO num2str(results.Pval(1,nn),'%10.3f') ' +/- ' num2str(results.PE(nn),'%10.3f')]);
    end
    fprintf(fid,'\n');
    
    if data.SHarmonics
        fprintf(fid,['Forcing Frequency (Hz):' pp.delimO num2str(data.FH) '\n']);
        fprintf(fid,['Forcing Strouhal Number:' pp.delimO num2str(data.StDF) '\n']);

        if pp.DT~=0
            fprintf(fid,'OASPL_Detoned (dB):'); fprintf(fid,'\t%.3f',results.OASPL_DT); fprintf(fid,'\n');
            fprintf(fid,['Average Energy_Detoned (dB):' pp.delimO num2str(results.AAE_DT,'%.3f') '\n']);
            fprintf(fid,['Normalized Jet Power_Detoned (dB):' pp.delimO num2str(results.NJP_DT,'%.3f') '\n']);     
        end
        
        if pp.nBaselines ~= 0 
            fprintf(fid,'dOASPL (dB):'); fprintf(fid,'\t%.4f',results.dOASPL); fprintf(fid,'\n');
            fprintf(fid,['dAAE (dB):' pp.delimO num2str(results.dAAE,'%.4f') '\n']);
            fprintf(fid,['dNJP (dB):' pp.delimO num2str(results.dNJP,'%.4f') '\n']);
            
            if pp.DT~=0
                fprintf(fid,'dOASPL_Detoned (dB):'); fprintf(fid,'\t%.4f',results.dOASPL_DT); fprintf(fid,'\n');
                fprintf(fid,['dAAE_Detoned (dB):' pp.delimO num2str(results.dAAE_DT,'%.4f') '\n']);
                fprintf(fid,['dNJP_Detoned (dB):' pp.delimO num2str(results.dNJP_DT,'%.4f') '\n']);
            end
        end
    else
        fprintf(fid,'OASPL_Smoothed (dB):'); fprintf(fid,'\t%.3f',results.OASPL_DT); fprintf(fid,'\n');
        fprintf(fid,['Average Energy_Smoothed (dB):' pp.delimO num2str(results.AAE_DT,'%.3f') '\n']);
        fprintf(fid,['Normalized Jet Power_Smoothed (dB):' pp.delimO num2str(results.NJP_DT,'%.3f') '\n']);
    end
    fclose(fid);
end