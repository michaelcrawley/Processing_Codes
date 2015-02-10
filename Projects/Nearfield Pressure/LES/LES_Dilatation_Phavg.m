function LES_Dilatation_Phavg(src)

    %Grab grid information
    grd = getfiles('grid.mat',src);
    load([src,filesep,grd{1}]);
    
    %Grab time and phase-averaged files
    flist = getfiles('mnf*.mat',src);
    for n = 1:length(flist)
        load([src,filesep,flist{n}],'u','v','w','u_r','u_theta');
        
        N = size(u);
        cart.C2 = zeros(N);
        cyl.C2 = zeros(N);
        for q = 1:N(3)
            cart.C2(:,:,q) = divergence(grid.y(:,:,q),grid.x(:,:,q),v(:,:,q),u(:,:,q)); %I have to switch the order (y,x instead of x,y) because of the direction that Matlab computes the derivatives in
            cyl.C2(:,:,q) = divergence(grid.r(:,:,q),grid.x(:,:,q),u_r(:,:,q),u(:,:,q)) + u_r(:,:,q)./grid.r(:,:,q); %divergence just calculates in the cartesian coordinate system
        end
        cyl.C3 = divergence(grid.r,grid.x,grid.theta,u_r,u,u_theta./grid.r) + u_r./grid.r;
        
        save([src,filesep,'DivCart',flist{n}],'-struct','cart');
        save([src,filesep,'DivCyl',flist{n}],'-struct','cyl');
    end    
end