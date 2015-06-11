function GenerateFigureFile(filename,data,pp,fig_dir)  
%Generates *.fig and *.png files of the spectra, which are saved in the
%directory specified by fig_dir.

%Last updated by Michael Crawley on 2011-09-08

    figure  %Plots all spectra on a single figure
	CM = colormap('jet');
    semilogx(data.Std,data.dBCorrected(:,1),'Color',CM(round(64/pp.Nch*1),:));
    LBL = 'Ch 1';  %Variable contains labels which will become column headers on output files
    hold on
    for m = 2:pp.Nch
        semilogx(data.Std,data.dBCorrected(:,m),'Color',CM(round(64/pp.Nch*m),:));
        LBL = [LBL '\tCh ' num2str(m)]; %#ok<AGROW>
    end
    LEGN = mat2cell(pp.ChPol,1,ones(1,pp.Nch));
    for m = 1:pp.Nch   %Plots smoothed spectra on top of original spectra
        semilogx(data.Std,data.dBCorrected_DT(:,m),'--k');
        LEGN{m} = num2str(LEGN{m});
    end
    legend(LEGN,'Location','NorthWest');
    xlabel(['Strouhal Number for D = ' num2str(pp.D) ' m Jet at M = ' num2str(pp.M)]);
    ylabel('SPL (dB)');
    title({['Average Spectrum for: ' filename],['Data Normalized for ' num2str(pp.NormD) ' x/D'],['Reynolds Number: ' num2str(data.ReN,'%.0f')]},'Interpreter','none');
    grid on;
    if length(pp.fLimits)==2
        xlim(pp.fLimits);
    else
        axis(pp.fLimits);
    end
    saveas(gcf,[fig_dir '\' filename(1:end-4),'.fig']); 
%     saveas(gcf,[fig_dir '\' filename(1:end-4),'.png']); 
%     print(gcf,'-dpng','-r300',[fig_dir '\' filename(1:end-4),'.png']);
    saveFigure_v2(gcf,[fig_dir '\' filename(1:end-4)],300)  
    close    %Closes figure
end