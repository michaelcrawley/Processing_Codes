function [] = Project3(dx,dt,Mx)
    %Euler equation solver, completed for AAE 694 Project 3 by 
    %Michael Crawley
    %Inputs:
    %       dx: spatial step size
    %       dt: temporal step size
    %       Mx: mean flow in x-Dfrection   
    tic;
    dy = dx;
    Nghost = 1;
    
    %Set up grid
    X = -100:dx:100;
    Y = 0:dy:100; %add ghost point to domain
    [x, y] = meshgrid(X,Y);
    foldername = [num2str(length(X)),' x ',num2str(length(Y)), ' Mx', strrep(num2str(Mx),'.','_')];
    [Ny, Nx] = size(x);
    
    %Set maximum time for simulation to run
    plot_times = [5 15 25 40 50 60 70 80];
    time_itr_max = max(plot_times);
    imax = length(0:dt:time_itr_max)-1;
    
    %Boundary conDftion coefficients
    r = sqrt(x.^2+y.^2);
    V = Mx*(x./r)+sqrt(1-(Mx*(y./r)).^2);
    A = V.*(x./r);
    B = V.*(y./r);
    C = V./(2*r); 
    
    %3rd order optimized time Dfscretization coefficients (b0, b1, b2, b3)
    b = [2.3025580883830 -2.4910075998482 1.5743409331815 -0.38589142217162];

    %initialize variables
    rho = zeros(length(Y),length(X),4);
    u = rho; v = rho; p = zeros(length(Y)+Nghost,length(X),4); %add ghost point to pressure domain    
    Drho = rho;
    Du = u;
    Dv = v;
    Dp = p;
    pswitch = p(:,:,1);
    
    %initial conditions
    rho(:,:,4) = 0.01*exp(-log(2)*(x.^2+(y-20).^2)/9);
    p(2:end,:,4) = 0.01*exp(-log(2)*(x.^2+(y-20).^2)/9);
    coefs = DRPlookup(1,5,dy)'; %create DRP coefficients to compute ghost value
    p(1,:,4) = (-(1/coefs(1))*coefs(2:end)*p(2:7,:,4));
    
    %create numerical differentiation matrices and cells
    mDx = mNumericalDerivative(1,6,dx,Nx,'-DRP');
    mDy = mNumericalDerivative(1,6,dy,Ny,'-DRP')';
    mDyp = mNumericalDerivative(1,6,dx,Ny+Nghost,'-DRP')';
    dqdx = cell(4,1); dqdy = dqdx;
    
    if exist(foldername,'file') ~= 7
        mkdir(pwd,foldername);
    end
    cd(foldername);
    for i = 1:imax %time step loop
        %time shift spatial derivative matrices
        Drho(:,:,1:3) = Drho(:,:,2:4);
        Du(:,:,1:3) = Du(:,:,2:4);
        Dv(:,:,1:3) = Dv(:,:,2:4);
        Dp(:,:,1:3) = Dp(:,:,2:4);
        
        dqdx{1} = rho(:,:,4)*mDx;
        dqdy{1} = mDy*rho(:,:,4);
        dqdx{2} = u(:,:,4)*mDx;
        dqdy{2} = mDy*u(:,:,4);
        dqdx{3} = v(:,:,4)*mDx;
        dqdy{3} = mDy*v(:,:,4);
        dqdx{4} = p(:,:,4)*mDx;
        dqdy{4} = mDyp*p(:,:,4);
        
        %interior nodes
        Drho(:,:,4) = -(Mx*dqdx{1}+dqdx{2})-dqdy{3}; 
        Du(:,:,4) = -(Mx*dqdx{2}+dqdx{4}(2:end,:));
        Dv(:,:,4) = -(Mx*dqdx{3})-dqdy{4}(2:end,:);
        Dp(2:end,:,4) = -(Mx*dqdx{4}(2:end,:)+dqdx{2})-dqdy{3};
        
        %left boundary (radiation)
        Drho(:,1:3,4) = A(:,1:3).*(-dqdx{1}(:,1:3))+B(:,1:3).*(-dqdy{1}(:,1:3))-C(:,1:3).*rho(:,1:3,4); 
        Du(:,1:3,4) = A(:,1:3).*(-dqdx{2}(:,1:3))+B(:,1:3).*(-dqdy{2}(:,1:3))-C(:,1:3).*u(:,1:3,4);
        Dv(:,1:3,4) = A(:,1:3).*(-dqdx{3}(:,1:3))+B(:,1:3).*(-dqdy{3}(:,1:3))-C(:,1:3).*v(:,1:3,4);
        Dp(2:end,1:3,4) = A(:,1:3).*(-dqdx{4}(2:end,1:3))+B(:,1:3).*(-dqdy{4}(2:end,1:3))-C(:,1:3).*p(2:end,1:3,4);
        
        %right boundary (outflow)
        Dp(2:end,end-2:end,4) = A(:,end-2:end).*(-dqdx{4}(2:end,end-2:end))+B(:,end-2:end).*(-dqdy{4}(2:end,end-2:end))-C(:,end-2:end).*p(2:end,end-2:end,4);
        Drho(:,end-2:end,4) = Mx*(-dqdx{1}(:,end-2:end)+dqdx{4}(2:end,end-2:end))+A(:,end-2:end).*(-dqdx{4}(2:end,end-2:end))+B(:,end-2:end).*(-dqdy{4}(2:end,end-2:end))-C(:,end-2:end).*p(2:end,end-2:end,4);
        
        %top boundary (radiation)        
        Drho(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{1}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{1}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*rho(end-2:end,4:end-3,4);
        Du(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{2}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{2}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*u(end-2:end,4:end-3,4);
        Dv(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{3}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{3}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*v(end-2:end,4:end-3,4);
        Dp(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{4}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{4}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*p(end-2:end,4:end-3,4);
               
        %calculate physical variables at new time period
        rhoswitch = rho(:,:,4)+dt*(b(1)*Drho(:,:,4)+b(2)*Drho(:,:,3)+b(3)*Drho(:,:,2)+b(4)*Drho(:,:,1));
        uswitch = u(:,:,4)+dt*(b(1)*Du(:,:,4)+b(2)*Du(:,:,3)+b(3)*Du(:,:,2)+b(4)*Du(:,:,1));
        vswitch = v(:,:,4)+dt*(b(1)*Dv(:,:,4)+b(2)*Dv(:,:,3)+b(3)*Dv(:,:,2)+b(4)*Dv(:,:,1));
        pswitch(2:end,:) = p(2:end,:,4)+dt*(b(1)*Dp(2:end,:,4)+b(2)*Dp(2:end,:,3)+b(3)*Dp(2:end,:,2)+b(4)*Dp(2:end,:,1));
        pswitch(1,:) = (-(1/coefs(1))*coefs(2:end)*p(2:7,:,4));
   
        %%shift physical variable matrices
        rho(:,:,1:3) = rho(:,:,2:4); 
        u(:,:,1:3) = u(:,:,2:4);
        v(:,:,1:3) = v(:,:,2:4);
        p(:,:,1:3) = p(:,:,2:4);
        
        %update physical variables
        rho(:,:,4) = rhoswitch;        
        u(:,:,4) = uswitch;        
        v(:,:,4) = vswitch;        
        p(:,:,4) = pswitch;        
        
        %plot variables at specified times
        if any(abs(dt*i - plot_times) <= eps) || i ==1
            pcolor(X,Y,rho(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Density contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            saveas(gcf,['density_i',num2str(i)],'fig');saveas(gcf,['density_i',num2str(i)],'png');
            
            pcolor(X,Y,u(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['U velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            saveas(gcf,['u_i',num2str(i)],'fig');saveas(gcf,['u_i',num2str(i)],'png');
            
            pcolor(X,Y,v(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['V velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            saveas(gcf,['v_i',num2str(i)],'fig');saveas(gcf,['v_i',num2str(i)],'png');
            
            pcolor(X,Y,p(2:end,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Pressure contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            saveas(gcf,['p_i',num2str(i)],'fig');saveas(gcf,['p_i',num2str(i)],'png');
            close all;
        end
    end
    cd ..
    compute_time = toc
end