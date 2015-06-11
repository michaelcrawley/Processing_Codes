function [] = Project2R2(dx,dt,Mx,My)
    %Euler equation solver, completed for AAE 694 Project 2 by 
    %Michael Crawley
    %Inputs:
    %       dx: spatial step size
    %       dt: temporal step size
    %       Mx: mean flow in x-Dfrection
    %       My: mean flow in y-Dfrection    
    tic;
    dy = dx;
    
    %Set up grid
    X = -100:dx:100;
    Y = -100:dy:100;
    [x, y] = meshgrid(X,Y);
    foldername = [num2str(length(X)),' x ',num2str(length(Y)), ' Mx', strrep(num2str(Mx),'.','_'), ' My',strrep(num2str(My),'.','_')];
    [~, N] = size(x);
    
    %Set maximum time for simulation to run
    time_itr_max = 0.05*12000;
    plot_times = 0.05*[1200 1600 2000 12000];
    imax = length(0:dt:time_itr_max)-1;
    
    %Boundary conDftion coefficients
    r = sqrt(x.^2+y.^2);
    V = Mx*(x./r)+My*(y./r)+sqrt(1-(Mx*(y./r)-My*(x./r)).^2);
    A = V.*(x./r);
    B = V.*(y./r);
    C = V./(2*r); 
    
    %3rd order optimized time Dfscretization coefficients (b0, b1, b2, b3)
    b = [2.3025580883830 -2.4910075998482 1.5743409331815 -0.38589142217162];

    %initialize variables
    rho = zeros(length(y),length(x),4);
    u = rho; v = rho; p = rho;    
    Drho = rho;
    Du = u;
    Dv = v;
    Dp = p;
    
    %initial conditions
    rho(:,:,4) = 0.01*exp(-log(2)*(x.^2+y.^2)/9)+0.001*exp(-log(2)*((x-67).^2+y.^2)/25);
    u(:,:,4) = 0.0004*y.*exp(-log(2)*((x-67).^2+y.^2)/25);
    v(:,:,4) = -0.0004*(x-67).*exp(-log(2)*((x-67).^2+y.^2)/25);
    p(:,:,4) = 0.01*exp(-log(2)*(x.^2+y.^2)/9);
    
    %create numerical differentiation matrices
    Mf = mNumericalDerivative(1,6,dx,N,'-DRP');
    Mft = Mf';
    dqdx = cell(4,1); dqdy = dqdx;
    
    for i = 1:imax %time step loop
        %time shift spatial derivative matrices
        Drho(:,:,1:3) = Drho(:,:,2:4);
        Du(:,:,1:3) = Du(:,:,2:4);
        Dv(:,:,1:3) = Dv(:,:,2:4);
        Dp(:,:,1:3) = Dp(:,:,2:4);
        
        dqdx{1} = rho(:,:,4)*Mf;
        dqdy{1} = Mft*rho(:,:,4);
        dqdx{2} = u(:,:,4)*Mf;
        dqdy{2} = Mft*u(:,:,4);
        dqdx{3} = v(:,:,4)*Mf;
        dqdy{3} =  Mft*v(:,:,4);
        dqdx{4} = p(:,:,4)*Mf;
        dqdy{4} =  Mft*p(:,:,4);
        
        %interior nodes
        Drho(:,:,4) = -(Mx*dqdx{1}+dqdx{2})-(My*dqdy{1}+dqdy{3}); 
        Du(:,:,4) = -(Mx*dqdx{2}+dqdx{4})-(My*dqdy{2});
        Dv(:,:,4) = -(Mx*dqdx{3})-(My*dqdy{3}+dqdy{4});
        Dp(:,:,4) = -(Mx*dqdx{4}+dqdx{2})-(My*dqdy{4}+dqdy{3});
        
%         %left boundary (radiation)
%         Drho(:,1:3,4) = A(:,1:3).*(-dqdx{1}(:,1:3))+B(:,1:3).*(-dqdy{1}(:,1:3))-C(:,1:3).*rho(:,1:3,4); 
%         Du(:,1:3,4) = A(:,1:3).*(-dqdx{2}(:,1:3))+B(:,1:3).*(-dqdy{2}(:,1:3))-C(:,1:3).*u(:,1:3,4);
%         Dv(:,1:3,4) = A(:,1:3).*(-dqdx{3}(:,1:3))+B(:,1:3).*(-dqdy{3}(:,1:3))-C(:,1:3).*v(:,1:3,4);
%         Dp(:,1:3,4) = A(:,1:3).*(-dqdx{4}(:,1:3))+B(:,1:3).*(-dqdy{4}(:,1:3))-C(:,1:3).*p(:,1:3,4);
%         
%         %right boundary
%         Dp(:,end-2:end,4) = A(:,end-2:end).*(-dqdx{4}(:,end-2:end))+B(:,end-2:end).*(-dqdy{4}(:,end-2:end))-C(:,end-2:end).*p(:,end-2:end,4);
%         if Mx == 0 %radiation boundary conditions            
%             Drho(:,end-2:end,4) = A(:,end-2:end).*(-dqdx{1}(:,end-2:end))+B(:,end-2:end).*(-dqdy{1}(:,end-2:end))-C(:,end-2:end).*rho(:,end-2:end,4);           
%             Du(:,end-2:end,4) = A(:,end-2:end).*(-dqdx{2}(:,end-2:end))+B(:,end-2:end).*(-dqdy{2}(:,end-2:end))-C(:,end-2:end).*u(:,end-2:end,4);
%             Dv(:,end-2:end,4) = A(:,end-2:end).*(-dqdx{3}(:,end-2:end))+B(:,end-2:end).*(-dqdy{3}(:,end-2:end))-C(:,end-2:end).*v(:,end-2:end,4); 
%         else %outflow boundary conditions
%             Drho(:,end-2:end,4) = Mx*(-dqdx{1}(:,end-2:end)+dqdx{4}(:,end-2:end))+A(:,end-2:end).*(-dqdx{4}(:,end-2:end))+My*(-dqdy{1}(:,end-2:end)+dqdy{4}(:,end-2:end))+B(:,end-2:end).*(-dqdy{4}(:,end-2:end))-C(:,end-2:end).*p(:,end-2:end,4);
%         end  
%         
%         %bottom boundary (radiation)
%         Drho(1:3,4:end-3,4) = A(1:3,4:end-3).*(-dqdx{1}(1:3,4:end-3))+B(1:3,4:end-3).*(-dqdy{1}(1:3,4:end-3))-C(1:3,4:end-3).*rho(1:3,4:end-3,4);      
%         Du(1:3,4:end-3,4) = A(1:3,4:end-3).*(-dqdx{2}(1:3,4:end-3))+B(1:3,4:end-3).*(-dqdy{2}(1:3,4:end-3))-C(1:3,4:end-3).*u(1:3,4:end-3,4);
%         Dv(1:3,4:end-3,4) = A(1:3,4:end-3).*(-dqdx{3}(1:3,4:end-3))+B(1:3,4:end-3).*(-dqdy{3}(1:3,4:end-3))-C(1:3,4:end-3).*v(1:3,4:end-3,4);
%         Dp(1:3,4:end-3,4) = A(1:3,4:end-3).*(-dqdx{4}(1:3,4:end-3))+B(1:3,4:end-3).*(-dqdy{4}(1:3,4:end-3))-C(1:3,4:end-3).*p(1:3,4:end-3,4);
%         
%         %top boundary
%         Dp(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{4}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{4}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*p(end-2:end,4:end-3,4);
%         if My == 0 %radiation boundary conditions
%             Drho(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{1}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{1}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*rho(end-2:end,4:end-3,4);
%             Du(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{2}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{2}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*u(end-2:end,4:end-3,4);
%             Dv(end-2:end,4:end-3,4) = A(end-2:end,4:end-3).*(-dqdx{3}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{3}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*v(end-2:end,4:end-3,4);
%         else %outflow boundary conditions          
%             Drho(end-2:end,4:end-3,4) = Mx*(-dqdx{1}(end-2:end,4:end-3)+dqdx{4}(end-2:end,4:end-3))+A(end-2:end,4:end-3).*(-dqdx{4}(end-2:end,4:end-3))+My*(-dqdy{1}(end-2:end,4:end-3)+dqdy{4}(end-2:end,4:end-3))+B(end-2:end,4:end-3).*(-dqdy{4}(end-2:end,4:end-3))-C(end-2:end,4:end-3).*p(end-2:end,4:end-3,4);
%         end       
        
        %calculate physical variables at new time period
        rhoswitch = rho(:,:,4)+dt*(b(1)*Drho(:,:,4)+b(2)*Drho(:,:,3)+b(3)*Drho(:,:,2)+b(4)*Drho(:,:,1));
        uswitch = u(:,:,4)+dt*(b(1)*Du(:,:,4)+b(2)*Du(:,:,3)+b(3)*Du(:,:,2)+b(4)*Du(:,:,1));
        vswitch = v(:,:,4)+dt*(b(1)*Dv(:,:,4)+b(2)*Dv(:,:,3)+b(3)*Dv(:,:,2)+b(4)*Dv(:,:,1));
        pswitch = p(:,:,4)+dt*(b(1)*Dp(:,:,4)+b(2)*Dp(:,:,3)+b(3)*Dp(:,:,2)+b(4)*Dp(:,:,1));                                        
   
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
        
%         plot variables at specified times
        if any(abs(dt*i - plot_times) <= eps) || i ==1
            figure;pcolor(X,Y,rho(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Density contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            
            figure;pcolor(X,Y,u(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['U velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            
            figure;pcolor(X,Y,v(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['V velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            
            figure;pcolor(X,Y,p(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Pressure contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
        end
    end
    compute_time = toc
end