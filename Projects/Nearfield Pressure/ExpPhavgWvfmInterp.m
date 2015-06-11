function ExpPhavgWvfmInterp(fname,DesGrid,save_name)


    %Load in phase-averaged waveforms from experimental data
    load(fname);
    
    %Get size info
    Lr = length(pf.sm_wvfm);
    Lt = min(cellfun(@length,pf.sm_wvfm));
    Lx = length(pf.pp.NFCh);
    Lff = length(pf.pp.FFCh);
    
    %Initialize Containers
    wvfm = zeros(Lr,Lx,Lt);
    ff_wvfm = zeros(Lr,Lff,Lt);
    X = zeros(Lr,Lx);
    Y = X;
    
    %Pull Phase averaged waveform and grid info
    FS = pf.pp.FS;
    for n = 1:Lr
        wvfm(n,:,:) = permute(pf.sm_wvfm{n}(1:Lt,pf.pp.NFCh),[3 2 1]);
        ff_wvfm(n,:,:) = permute(pf.sm_wvfm{n}(1:Lt,pf.pp.FFCh),[3 2 1]);
        X(n,:) = pf.phys{n}.x;
        Y(n,:) = pf.phys{n}.y;
    end
    ff_wvfm = squeeze(mean(ff_wvfm,1));
    
    %Reshape Grid, if necessary
    [X,Y,wvfm] = ReshapeGrid(X,Y,wvfm);
    
    %Interpolate onto desired grid
    for n = 1:Lt
%         vals = wvfm(:,:,n);
%         F = TriScatteredInterp(X(:),Y(:),vals(:));
        output.p(:,:,n) = griddata(X,Y,wvfm(:,:,n),DesGrid.x,DesGrid.y,'linear');
%         output.p(:,:,n) = F(DesGrid.x,DesGrid.y);
    end
    
    %Get misc info
    output.x = DesGrid.x;
    output.y = DesGrid.y;
    output.dt = 1/FS;
    output.ffp = ff_wvfm;
    
    %Save mat file
    [~,~,ext] = fileparts(save_name);
    if strcmpi(ext,''), save_name = [save_name '.mat']; end
    save(save_name,'-struct','output');
end