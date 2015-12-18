function AcousticPlotter(src_dir,flist,varargin)
    cmds = varargin;

    if ~exist('varargin','var')
       error('Error: no forcing options were chosen'); 
    end
    
    %%plot individual spectra versus baselines
    if any(strcmpi(varargin,'-Spectra'))
        bflist = varargin{1};
        Spectra(src_dir,flist,bflist);
    end
    
    %%Plots individual baseline spectra
    if any(strcmpi(varargin,'-bSpectra'))
        bSpectra(src_dir,flist);        
    end
    
    %%Plots all baseline spectra
    if any(strcmpi(varargin,'-abSpectra'))
        abSpectra(src_dir,flist);        
    end
    
    %%plots delta OASPL
    if any(strncmpi(varargin,'-doaspl',7))
        doaspl(src_dir,flist,cmds);
    end
    
    %%Create plot of baseline OASPL values - run number versus OASPL for all polar angle
    if any(strncmpi(varargin,'-boaspl',6))
        boaspl(src_dir,flist,varargin);
    end
    
    %%plot Nonlinear metrics for forced cases
    if strcmpi(varargin{1},'-NLM')
        NLM(src_dir,flist,varargin);
    end
    
    %%plot Nonlinear metrics for baseline cases
    if strcmpi(varargin{1},'-bNLM')
        bNLM(src_dir,flist,varargin);
    end
    
    %%plots delta OASPL standard deviation (only for averaged mat files)
    if any(strncmpi(varargin,'-doasplstd',10))
        doasplstd(src_dir,flist,varargin);
    end
    
    %%plot Nonlinear metrics standard deviation for forced cases (only for
    %%averaged mat files)
    if strcmpi(varargin{1},'-NLMstd')
        NLMstd(src_dir,flist,varargin);
    end
end

function NLMstd(src_dir, flist, varargin)
    pull = varargin{2};
    Stdf = zeros(1,length(flist));
    signal = zeros(length(flist),12); %fix for variable number of microphones

    %pull data from mat files
    for n = 1:length(flist) 
       tmp = load(fullfile(src_dir,flist{n}));
       Stdf(n) = eval(['tmp.',flist{n}(1:end-4),'.data.StDF']);
       signal(n,:) = eval(['tmp.',flist{n}(1:end-4),'.results.', pull,'std']);
    end
    ChPol = eval(['tmp.',flist{n}(1:end-4),'.processing_params.ChPol']);

    %create plot        
    if length(varargin) == 2
        filename = input('Input name for contour plot: ','s');
        plottitle = '';
    elseif length(varargin) == 3
        filename = varargin{3};
        plottitle = '';
    elseif length(varargin) == 4
        filename = varargin{3};
        plottitle = varargin{4};
    end

    H = figure;
    pcolor(ChPol,Stdf,signal); shading interp; colorbar;
    xlabel('Polar Angle (degrees)'); ylabel('St_D_F'); title(plottitle,'Interpreter','none');
    set(gca,'XDir','reverse');
    set(gca,'Xtick',30:10:130);

    %save and close plot
    saveas(H,filename,'fig'); 
    saveFigure_v2(H,filename,600); 
    close(H);
end

function doasplstd(src_dir, flist, varargin)
    Stdf = zeros(1,length(flist));
    doasplstd = zeros(length(flist),12); %fix for variable number of microphones

    %pull data from mat files
    for n = 1:length(flist) 
       tmp = load(fullfile(src_dir,flist{n}));
       Stdf(n) = eval(['tmp.',flist{n}(1:end-4),'.data.StDF']);
       doasplstd(n,:) = eval(['tmp.',flist{n}(1:end-4),'.results.dOASPL_DTstd']);
    end
    ChPol = eval(['tmp.',flist{n}(1:end-4),'.processing_params.ChPol']);

    %create plot
    cmdstr = varargin{strncmpi(varargin,'-doasplstd',10)};
    I = regexp(cmdstr,'-');

    if length(I) == 1
        filename = input('Input name for contour plot: ','s');
        plottitle = '';
    elseif length(I) == 2
        filename = cmdstr(I(2)+1:end);
        plottitle = '';
    else
        filename = cmdstr(I(2)+1:I(3)-1);
        plottitle = cmdstr(I(3)+1:end);
    end

    H = figure;
    [cs,h]=contourf(ChPol,Stdf,doasplstd); colorbar; clabel(cs,h,'rotation',0); 
    xlabel('Polar Angle (degrees)'); ylabel('St_D_F'); title(plottitle,'Interpreter','none');
    set(gca,'XDir','reverse');
    set(gca,'Xtick',30:10:130);        
    set(gca,'XTickmode','manual');
    grid('on');

    %save and close plot
    saveas(H,filename,'fig'); 
    saveFigure_v2(H,filename,600); 
    close(H);
end

function bNLM(src_dir,flist,varargin)
    pull = varargin{2};
    tnum = zeros(1,length(flist));
    signal = zeros(length(flist),12); %fix for variable number of microphones

    %pull data from mat files
    for n = 1:length(flist) 
        tmp = load(fullfile(src_dir,flist{n}));
        tnum(n) = eval(['tmp.',flist{n}(1:end-4),'.data.tnum']);
        signal(n,:) = eval(['tmp.',flist{n}(1:end-4),'.data.', pull]);
    end
    [~,I] = sort(tnum);
    signal = signal(I,:);
    ChPol = eval(['tmp.',flist{n}(1:end-4),'.processing_params.ChPol']);

    %create plot        
    if length(varargin) == 2
        filename = input('Input name for contour plot: ','s');
        plottitle = '';
    elseif length(varargin) == 3
        filename = varargin{3};
        plottitle = '';
    elseif length(varargin) == 4
        filename = varargin{3};
        plottitle = varargin{4};
    end

    H = figure;
    pcolor(ChPol,1:length(tnum),signal); shading interp; colorbar;
    xlabel('Polar Angle (degrees)'); ylabel('set'); title(plottitle,'Interpreter','none');
    set(gca,'XDir','reverse');
    set(gca,'Xtick',30:10:130);

    %save and close plot
    saveas(H,filename,'fig'); 
    saveFigure_v2(H,filename,600); 
        close(H);
end

function NLM(src_dir,flist,varargin)
    pull = varargin{2};
    Stdf = zeros(1,length(flist));
    signal = zeros(length(flist),12); %fix for variable number of microphones

    %pull data from mat files
    for n = 1:length(flist) 
       tmp = load(fullfile(src_dir,flist{n}));
       Stdf(n) = eval(['tmp.',flist{n}(1:end-4),'.data.StDF']);
       signal(n,:) = eval(['tmp.',flist{n}(1:end-4),'.results.', pull]);
    end
    ChPol = eval(['tmp.',flist{n}(1:end-4),'.processing_params.ChPol']);

    %create plot        
    if length(varargin) == 2
        filename = input('Input name for contour plot: ','s');
        plottitle = '';
    elseif length(varargin) == 3
        filename = varargin{3};
        plottitle = '';
    elseif length(varargin) == 4
        filename = varargin{3};
        plottitle = varargin{4};
    end

    H = figure;
    pcolor(ChPol,Stdf,signal); shading interp; colorbar;
    xlabel('Polar Angle (degrees)'); ylabel('St_D_F'); title(plottitle,'Interpreter','none');
    set(gca,'XDir','reverse');
    set(gca,'Xtick',30:10:130);

    %save and close plot
    saveas(H,filename,'fig'); 
    saveFigure_v2(H,filename,600); 
    close(H);
end

function boaspl(src_dir,flist,varargin)
    oaspl = zeros(length(flist),12); %fix for variable number of microphones

    %pull data from mat files
    for n = 1:length(flist) 
       tmp = load(fullfile(src_dir,flist{n}));
       time(n) = eval(['tmp.',flist{n}(1:end-4),'.data.tnum']);
       oaspl(n,:) = eval(['tmp.',flist{n}(1:end-4),'.results.OASPL']);
    end
    ChPol = eval(['tmp.',flist{n}(1:end-4),'.processing_params.ChPol']);
    [~,I] = sort(time);
    oaspl = oaspl(I,:); %resort in ascending order of time recorded

    %create plot
    cmdstr = varargin{strncmpi(varargin,'-boaspl',7)};
    I = regexp(cmdstr,'-');

    if length(I) == 1
        filename = input('Input name for contour plot: ','s');
        plottitle = '';
    elseif length(I) == 2
        filename = cmdstr(I(2)+1:end);
        plottitle = '';
    else
        filename = cmdstr(I(2)+1:I(3)-1);
        plottitle = cmdstr(I(3)+1:end);
    end

    H = figure;
    plot(1:n,oaspl); xlabel('baseline set'); ylabel('OASPL (dB)'); legend(num2str(ChPol'),'Location','NorthEastOutside');

    %save and close plot
    saveas(H,filename,'fig'); 
    saveas(H,filename,'png'); 
    close(H);
end

function doaspl(src_dir,flist,cmds)

    Stdf = zeros(1,length(flist));
    doaspl = zeros(length(flist),12); %fix for variable number of microphones

    %pull data from mat files
    for n = 1:length(flist) 
       tmp = load(fullfile(src_dir,flist{n}));
       Stdf(n) = eval(['tmp.',flist{n}(1:end-4),'.data.StDF']);
       doaspl(n,:) = eval(['tmp.',flist{n}(1:end-4),'.results.dOASPL_DT']);
    end
    ChPol = eval(['tmp.',flist{n}(1:end-4),'.processing_params.ChPol']);

    %create plot
    cmdstr = cmds{strncmpi(cmds,'-doaspl',7)};
    I = regexp(cmdstr,'-');

    if length(I) == 1
        filename = input('Input name for contour plot: ','s');
        plottitle = '';
    elseif length(I) == 2
        filename = cmdstr(I(2)+1:end);
        plottitle = '';
    else
        filename = cmdstr(I(2)+1:I(3)-1);
        plottitle = cmdstr(I(3)+1:end);
    end

    H = figure;
    increment = 0.25;
    minval = floor(min(min(doaspl/increment)))*increment;
    maxval = ceil(max(max(doaspl/increment)))*increment;
    [cs,h]=contourf(ChPol,Stdf,doaspl,minval:increment:maxval); colorbar; clabel(cs,h,'rotation',0); 
    xlabel('Polar Angle (degrees)'); ylabel('St_D_F'); title(plottitle,'Interpreter','none');
    set(gca,'XDir','reverse');
    set(gca,'Xtick',30:10:130);

    set(gca,'XTickmode','manual');
    grid('on');
    load('Z:\My Documents\Matlab\General Use\Colormaps\cool_range_custom.mat');
    [map] = CreateCustomMap([min(min(doaspl)) maxval-increment],trange,brange);
    colormap(map);

    %save and close plot
    saveas(H,filename,'fig'); 
    saveFigure_v2(H,filename,600); 
    close(H);
end

function bSpectra(src_dir,flist)
    for n = 1:length(flist)
        [~,filename] = fileparts(flist{n});
        mkdir([src_dir filesep filename]);

        tmp = load(fullfile(src_dir,flist{n}));
        spectra = eval(['tmp.',filename,'.data.dBCorrected']);
        Std = eval(['tmp.',filename,'.data.Std']);
        ChPol = eval(['tmp.',filename,'.processing_params.ChPol']);  

        for nn = 1:length(ChPol) 
            h = figure;
            semilogx(Std,spectra(:,nn)); xlabel('St_d'); ylabel('SPL (dB)'); title([filename,' Polar Angle: ', num2str(ChPol(nn))],'Interpreter','none');
            xlim([0.01 6]);
            saveas(h,[src_dir filesep filename filesep 'Pol ' num2str(ChPol(nn))],'fig');
            saveas(h,[src_dir filesep filename filesep 'Pol ' num2str(ChPol(nn))],'png');
            close(h);
        end

    end
end

function Spectra(src_dir,flist,bflist)
    %load baseline data
    bspectra = cell(1,length(bflist));
    bT = zeros(1,length(bflist));
    bspectra = zeros(length(bflist),4086,12); %need to fix for variables
    bStd = zeros(length(bflist),4086);

    for n = 1:length(bflist)
        tmp = load(fullfile(src_dir,bflist{n}));
        bspectra(n,:,:) = eval(['tmp.',bflist{n}(1:end-4),'.data.dBCorrected']);
        bT(n) = eval(['tmp.',bflist{n}(1:end-4),'.data.To']);
        bStd(n,:) = eval(['tmp.',bflist{n}(1:end-4),'.data.Std']);
    end
    [bT,I] = sort(bT);
    bspectra = bspectra(I,:,:);
    bStd = bStd(I,:);

    %removes non-unique spectra based on stagnation temperature
    [bT,I] = unique(bT);
    bspectra = bspectra(I,:,:);
    bStd = bStd(I,:);        

    for n = 1:length(flist)
        [~,filename] = fileparts(flist{n});
        mkdir([src_dir filesep filename]);

        tmp = load(fullfile(src_dir,flist{n}));
        spectra = eval(['tmp.',filename,'.data.dBCorrected']);
        Std = eval(['tmp.',filename,'.data.Std']);
        To = eval(['tmp.',filename,'.data.To']); 
        if To > max(bT), To = max(bT); end
        ChPol = eval(['tmp.',filename,'.processing_params.ChPol']);  
        baselinespectra = squeeze(interp1(bT,bspectra,To,'linear','extrap'));
        baselineStd = squeeze(interp1(bT,bStd,To,'linear','extrap'));

        for nn = 1:length(ChPol) 
            h = figure;
            semilogx(baselineStd,baselinespectra(:,nn),'k',Std,spectra(:,nn)); xlabel('St_d'); ylabel('SPL (dB)'); title([filename,' Polar Angle: ', num2str(ChPol(nn))],'Interpreter','none');
            xlim([0.01 6]);
            legend('Baseline','Forced','Location','Northwest');
            saveas(h,[src_dir filesep filename filesep 'Pol ' num2str(ChPol(nn))],'fig');
            saveas(h,[src_dir filesep filename filesep 'Pol ' num2str(ChPol(nn))],'png');
            close(h);
        end

    end
end

function abSpectra(src_dir,flist)
    mkdir([src_dir filesep 'abSpectra']);
    N = length(flist);
    
    for n = 1:N
        [~,filename] = fileparts(flist{n});
        tmp = load(fullfile(src_dir,flist{n}));
        spectra(:,:,n) = eval(['tmp.',filename,'.data.dBCorrected']);
        Std(:,n) = eval(['tmp.',filename,'.data.Std']);
        ChPol = eval(['tmp.',filename,'.processing_params.ChPol']);
        Tj(n) = eval(['tmp.',filename,'.data.To']);
        tnum(n) = eval(['tmp.',filename,'.data.tnum']);
    end
    
    %reorder based on timestamp
    [~,I] = sort(tnum);
    spectra = spectra(:,:,I);
    Std = Std(:,I);
    Tj = Tj(I);
    
    %add offset
    offset = 5;
    for n = 1:N
        spectra(:,:,n) = spectra(:,:,n) + offset*(n-1);
    end
    
    %get x limits
    xmin = min(min(Std));
    xmax = max(max(Std));
        
    for n = 1:length(ChPol)
        h = figure;  
        CM = colormap('jet');
        for nn = 1:N
            semilogx(Std(:,nn),spectra(:,n,nn),'Color',CM(round(64/N*nn),:));
            hold on;
        end
        
        ymin = min(min(spectra(:,n,:)));
        ymax = max(max(spectra(:,n,:)));
        
        xlabel('St_D');
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        title(['\theta = ',num2str(ChPol(n)),'^o']);
        set(gca,'yticklabel','');
        legend(num2str(Tj'),'Location','Northwestoutside');
        saveas(h,[src_dir filesep 'abSpectra' filesep num2str(ChPol(n))],'fig');
        print('-dpng','-r600',[src_dir filesep 'abSpectra' filesep num2str(ChPol(n))]);
        close(h);
    end
    
end