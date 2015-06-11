function GenerateFFTFile(filename,data,pp,fft_dir) 
 %Generates *.fft and *.S.fft files based on values in 'pp' and 'data' structures.  
 %*.fft and *.S.fft files are saved in the directory specified by fft_dir. 
 
 %Last updated by Michael Crawley on 2011-07-04

    LBL = 'Ch 1';  %Variable contains labels which will become column headers on output files
    for m = 2:pp.Nch
        LBL = [LBL '\tCh ' num2str(m)]; %#ok<AGROW>
    end    

    fid = fopen([fft_dir '\' filename(1:end-3) 'fft' pp.EXT],'w'); %Opens .fft* file to write column headers 
    fprintf(fid,['Frequency (Hz)\tStrouhal Number' pp.delimO LBL '\n']);
    fclose(fid);
    dlmwrite([fft_dir '\' filename(1:end-3) 'fft' pp.EXT],[pp.fx data.Std data.dBCorrected], 'delimiter', pp.delimO, '-append');    %writes spectral data to file 
    
    fid = fopen([fft_dir '\' filename(1:end-3) 'S.fft' pp.EXT],'w'); %Opens .S.fft* file to write column headers 
    fprintf(fid,['Frequency (Hz)\tStrouhal Number' pp.delimO LBL '\n']);
    fclose(fid); clear LBL;
    dlmwrite([fft_dir '\' filename(1:end-3) 'S.fft' pp.EXT],[pp.fx data.Std data.dBCorrected_DT], 'delimiter', pp.delimO, '-append');  %writes smoothed spectral data to file 
end