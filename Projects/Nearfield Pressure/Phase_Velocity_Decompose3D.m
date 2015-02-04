function [subsonic supersonic PSD f intwvfm] = Phase_Velocity_Decompose3D(src_dir,flist,out_dir,a,theta)

    %Constants
    windowfun = @(x) tukeywin(x,0.1);
    Nx = 22; %number of axial points for upsampling
    D = 0.0254;
    FS = 2e5;
    dt = 1/FS;
    beta_const = 22/sqrt(-log(1e-6));
    NFCh = 4:19;
    transform = [secd(theta) -tand(theta)]; %transform vector to orthogonal component vectors
    
    for n = 1:length(flist)
        %Read in data and sort arrays
        data = load([src_dir filesep flist{n}]);
        Nt = min(cellfun(@(x) size(x,1),data.pf.sm_wvfm));
        sig = cellfun(@(x) x(1:Nt,:),data.pf.sm_wvfm,'uniformoutput',false);
        NFsig = cell2mat(cellfun(@(x) permute(x(:,NFCh),[3 2 1]),sig,'uniformoutput',false)');
        x = cell2mat(cellfun(@(X) X.x,data.pf.phys,'uniformoutput',false)')*D;
        y = cell2mat(cellfun(@(Y) Y.y,data.pf.phys,'uniformoutput',false)')*D;
        [X Y NFsig] = ReshapeGrid(x,y,NFsig);
        
        %Upsample for FFT
        Ny = size(X,1);
        intX = repmat(linspace(X(1),X(end),Nx),Ny,1);
        intY = repmat(Y(:,1),1,Nx)+(intX-intX(1))*tand(theta);
        intY(1,:) = intY(1,:)+5*eps;
        intY(end,:) = intY(end,:)-5*eps;
        intwvfm = zeros(Ny,Nx,Nt);
        for k = 1:Nt
            intwvfm(:,:,k) = griddata(X,Y,NFsig(:,:,k),intX,intY,'cubic');
        end      
        
        %Compute 3D FFT
        difx = diff(intX,1,2);
        dify = diff(intY,1,1);
        dx = mean(difx(:));
        dy = mean(dify(:));
        [PSD, f, S, xm] = PSDN(intwvfm,[1 2 3],[dy dx dt],windowfun,true);
        
        %Calculate Weigth Vectors
        ky = repmat(permute(f{1},[1 2 3]),[1 Nx Nt]);
        km = repmat(permute(f{2},[2 1 3]),[Ny 1 Nt]);
        omega = repmat(permute(f{3},[2 3 1]),[Ny Nx 1]);
        kx = transform(1)*km+transform(2)*ky;
        kb = sqrt(kx.^2+ky.^2); %assumes plane waves
        ka = omega/a;
        dky = ky(2);
        dkx = kx(2);
        dkax = abs(kx)-abs(ka);
        dkay = abs(ky)-abs(ka);
        betax = beta_const*dkx;
        betay = beta_const*dky;
        W = exp(-sqrt((dkax/betax).^2+(dkay/betay).^2));
        vel = abs(omega./kb); %phase velocity
        vel(isinf(vel) | isnan(vel)) = 0; %get rid of inf value at kb = 0
        supW = W;
        supW(vel > a) = 1;
        subW = ones(size(vel))-supW;        
        
        %Compute 3D iFFT
        subsonic = real(iPSDN(S.*subW,[1 2 3],xm));
        supersonic = real(iPSDN(S.*supW,[1 2 3],0));
        
        %Save Data (if desired)
        [~,fname] = fileparts(flist{n});
        if nargout == 0
            save([out_dir filesep fname ' 3D decomp.mat'],'a','windowfun','subsonic','supersonic','data','PSD','intwvfm','f','intX','intY','W','supW','subW');
        end
    end
end