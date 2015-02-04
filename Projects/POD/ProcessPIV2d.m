function ProcessPIV2d(src_dir,out_dir,To,M,D)

    str = genpath(src_dir);
    sd = regexp(str,';','split'); 
    test = regexp(sd,'PostProc', 'once');
    test2 = cellfun(@(cell) ~isempty(cell),test);
    subdir = sd(test2);
    l = length(src_dir);
    mkdir(out_dir);
    
    for n = 1:length(subdir)
        Tj = To(n)/(1+1/5*M^2);    %jet exit temp - K
       
        flist = getfiles('*.VC7',subdir{n});
        [flist X x y] = PreProcessPIV2d(subdir{n},flist); %#ok<ASGLU> %throw out garbage images
        x = x/1000/D; y = y/1000/D;
        [phi, lambda] = POD2d(X);
        P = length(x); N = length(y);
        U = reshape(X(1:P*N,:),P,N,[]);
        Um = mean(U,3);
        V = reshape(X(P*N+1:end,:),P,N,[]);
        Vm = mean(V,3); %#ok<NASGU>
        POD.U = reshape(phi(1:P*N,:),P,N,[]);
        POD.V = reshape(phi(P*N+1:end,:),P,N,[]);
        POD.eigs = lambda;
        
        FWHM = findFWHM(x,y,Um);    %#ok<NASGU>%jet width - y/D

        qx = find(x>1, 1 );
        qy = find(y>0, 1, 'last' );
        core = Um(1:qx,qy-1:qy+1);
        Uj = mean(core(:));  %jet exit velocity - m/s
        
        CL = Um(:,qy)/sqrt(1.4*287.05*Tj);  %#ok<NASGU>%jet centerline mach number - dimensionless
        TKE = 1/Uj^2*(std(U,0,3).^2 + std(V,0,3).^2);  %#ok<NASGU>%2D TKE - dimensionless
        
        %parse subdir name to find set name
        split = regexp(subdir{n}(l+2:end),filesep,'split');
        fname = split{1};
        save([out_dir filesep fname '.mat'],'flist','x','y','U','Um','V','Vm','POD','FWHM','Uj','CL','TKE');
    end

end