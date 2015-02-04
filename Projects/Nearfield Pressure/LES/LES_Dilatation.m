function LES_Dilatation(src)

    %Compute Dilatation in Cartesian coordinates
    ufile = getfiles('u.mat',src);
    u = load([src,filesep,ufile{1}]);
    grid = u.grid;
    t = u.t;
    trigger = u.trigger;
%     xcoefs = sparse(nuTSE(3,u.grid.x(1,:),1));
%     ycoefs = sparse(nuTSE(3,u.grid.y(:,1),1));
%     zcoefs = sparse(nuTSE(3,u.grid.z(:,1),1));  %can I even use this? should be numerical error....
%     rcoefs = sparse(nuTSE(3,u.grid.r(:,1),1));
%     thetacoefs = sparse(nuTSE(3,u.grid.theta(:,1),1));  %probably garbage

    %Process trigger in order to phase-average
    N = size(u.data);
    Nfp = unique(diff(find(diff(trigger) > 0))); %number of points in forcing period
    cut = mod(N(4),Nfp); %remove partial period at end of set
    nf = (N(4)-cut)/Nfp;
    
    cartDilatation = zeros(N); 
    PhavgcartDilatation = zeros([N(1:3),Nfp]);
    
    %partial u / partial x
    for n = 1:N(4)
        [px, ~] = gradient(u.data(:,:,1,n), grid.x(1,:), grid.y(:,1)); 
        cartDilatation(:,:,1,n) = cartDilatation(:,:,1,n) + px;
    end
    phavg_u = mean(reshape(u.data(:,:,1,1:end-cut),[N(1:3),Nfp,nf]),5);
    for n = 1:Nfp
        [px,~] = gradient(phavg_u(:,:,1,n),grid.x(1,:), grid.y(:,1)); 
        PhavgcartDilatation(:,:,1,n) = PhavgcartDilatation(:,:,1,n) + px;
    end    
    clear u;
    
    %partial v/ partial y
    vfile = getfiles('v.mat',src);
    v = load([src,filesep,vfile{1}],'data');
    for n = 1:N(4)
        [~,py] = gradient(v.data(:,:,1,n),grid.x(1,:),grid.y(:,1));
        cartDilatation(:,:,1,n) = cartDilatation(:,:,1,n) + py; 
    end
    phavg_v = mean(reshape(v.data(:,:,1,1:end-cut),[N(1:3),Nfp,nf]),5);
    for n = 1:Nfp
        [~,py] = gradient(phavg_v(:,:,1,n),grid.x(1,:), grid.y(:,1)); 
        PhavgcartDilatation(:,:,1,n) = PhavgcartDilatation(:,:,1,n) + py;
    end 
    clear v;
    
    %save out two-component dilatation...
    save([src,filesep,'cartDilatation2C.mat'],'cartDilatation','PhavgcartDilatation','grid','t','trigger');
    
    %partial w/ partial z
    wfile = getfiles('w.mat',src);
    w = load([src,filesep,wfile{1}],'data');
    for n = 1:N(4)
        [~,pz] = gradient(w.data(:,:,1,n),grid.x(1,:),grid.z(:,1));
        cartDilatation(:,:,1,n) = cartDilatation(:,:,1,n) + pz; %pre-multiply for column derivative....
    end
    phavg_w = mean(reshape(w.data(:,:,1,1:end-cut),[N(1:3),Nfp,nf]),5);
    for n = 1:Nfp
        [~,pz] = gradient(phavg_w(:,:,1,n),grid.x(1,:), grid.z(:,1)); 
        PhavgcartDilatation(:,:,1,n) = PhavgcartDilatation(:,:,1,n) + pz;
    end 
    clear w;
    
    %save out full dilatation field
    save([src,filesep,'cartDilatation3C.mat'],'cartDilatation','PhavgcartDilatation','grid','t','trigger');
    clear cartDilatation;
    
    %Compute Dilatation in Cylidrical coordinates
    cylDilatation = zeros(N);
    PhavgcylDilatation = zeros([N(1:3),Nfp]);
        
    %partial u / partial x
    ufile = getfiles('u.mat',src);
    u = load([src,filesep,ufile{1}],'data');
    for n = 1:N(4)
        [px, ~] = gradient(u.data(:,:,1,n), grid.x(1,:), grid.r(:,1)); 
        cylDilatation(:,:,1,n) = cylDilatation(:,:,1,n) + px;
    end
    for n = 1:Nfp
        [px,~] = gradient(phavg_u(:,:,1,n),grid.x(1,:), grid.y(:,1)); 
        PhavgcylDilatation(:,:,1,n) = PhavgcylDilatation(:,:,1,n) + px;
    end 
    clear u;
    
    %partial u_r / partial r
    u_rfile = getfiles('u_r.mat',src);
    u_r = load([src,filesep,u_rfile{1}],'data');
    for n = 1:N(4)
        [~,pr] = gradient(grid.r.*u_r.data(:,:,1,n),grid.x(1,:), grid.r(:,1));
        cylDilatation(:,:,1,n) = cylDilatation(:,:,1,n) + pr./grid.r; %pre-multiply for row derivative...
    end
    phavg_u_r = mean(reshape(u_r.data(:,:,1,1:end-cut),[N(1:3),Nfp,nf]),5);
    for n = 1:Nfp
        [~,py] = gradient(phavg_u_r(:,:,1,n),grid.x(1,:), grid.y(:,1)); 
        PhavgcartDilatation(:,:,1,n) = PhavgcartDilatation(:,:,1,n) + py;
    end
    clear u_r;
    
    %save out two-component dilatation...
    save([src,filesep,'cylDilatation2C.mat'],'cylDilatation','PhavgcylDilatation','grid','t','trigger');
    
    %partial u_theta / partial theta
    u_thetafile = getfiles('u_theta.mat',src);
    u_theta = load([src,filesep,u_thetafile{1}],'data');
    for n = 1:N(4)
        [~,ptheta] = gradient(u_theta.data(:,:,1,n)./grid.r,grid.x(1,:), grid.r(:,1));
        cylDilatation(:,:,1,n) = cylDilatation(:,:,1,n) + ptheta; %pre-multiply for row derivative...
    end
    phavg_u_theta = mean(reshape(u_theta.data(:,:,1,1:end-cut),[N(1:3),Nfp,nf]),5);
    for n = 1:Nfp
        [~,pz] = gradient(phavg_u_theta(:,:,1,n),grid.x(1,:), grid.z(:,1)); 
        PhavgcartDilatation(:,:,1,n) = PhavgcartDilatation(:,:,1,n) + pz;
    end
    clear u_theta;
    
    %save out full dilatation field
    save([src,filesep,'cylDilatation3C.mat'],'cylDilatation','PhavgcylDilatation','grid','t','trigger'); 
end