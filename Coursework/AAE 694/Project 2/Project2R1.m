function [] = Project2R1(dx,dt,Mx,My)
    %Euler equation solver, completed for AAE 694 Project 2 by 
    %Michael Crawley
    %Not complete; error in boundary
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
    [M, N] = size(x);   
    
    %Set maximum time for simulation to run
    time_itr_max = 0.05*5000;
    plot_times = 0.05*[1200 1600 2000 12000];
    imax = length(0:dt:time_itr_max)-1;
    
    %Init fields
    iField = {floor(M/4):ceil(3*M/4),floor(N/4):N};
    iSig = {1:M,10};
    Field = zeros(length(iField{1}),length(iField{2}),imax);
    Sig = zeros(M,imax);
    locations.Field.Y = Y(iField{1});
    locations.Field.X = X(iField{2});
    locations.Sig.Y = Y(iSig{1});
    locations.Sig.X = X(iSig{2});
    
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
    Mf = mNumericalDerivative(1,6,dx,N,'-DRP')';
    Mp = mNumericalDerivative(1,6,dx,7,'-DRP')';
    

    for i = 1:imax %time step loop
        %time shift spatial derivative matrices
        Drho(:,:,1:3) = Drho(:,:,2:4);
        Du(:,:,1:3) = Du(:,:,2:4);
        Dv(:,:,1:3) = Dv(:,:,2:4);
        Dp(:,:,1:3) = Dp(:,:,2:4);
        
        %interior nodes
        Drhoi = -(Mx*rho(:,:,4)+u(:,:,4))*Mf+(-(My*rho(:,:,4)+v(:,:,4))'*Mf)'; 
        Dui = -(Mx*u(:,:,4)+p(:,:,4))*Mf+(-(My*u(:,:,4))'*Mf)';
        Dvi = -(Mx*v(:,:,4))*Mf+(-(My*v(:,:,4)+p(:,:,4))'*Mf)';
        Dpi = -(Mx*p(:,:,4)+u(:,:,4))*Mf+(-(My*p(:,:,4)+v(:,:,4))'*Mf)';
        
        %left boundary (radiation)
        Drhol = A(:,1:7).*(-(rho(:,1:7,4))*Mp)+B(:,1:7).*(-(rho(:,1:7,4))'*Mf)'-C(:,1:7).*rho(:,1:7,4); 
        Dul = A(:,1:7).*(-(u(:,1:7,4))*Mp)+B(:,1:7).*(-(u(:,1:7,4))'*Mf)'-C(:,1:7).*u(:,1:7,4);
        Dvl = A(:,1:7).*(-(v(:,1:7,4))*Mp)+B(:,1:7).*(-(v(:,1:7,4))'*Mf)'-C(:,1:7).*v(:,1:7,4);
        Dpl = A(:,1:7).*(-(p(:,1:7,4))*Mp)+B(:,1:7).*(-(p(:,1:7,4))'*Mf)'-C(:,1:7).*p(:,1:7,4);
        
        %right boundary
        if Mx == 0 %radiation boundary conditions            
            Drhor = A(:,end-6:end).*(-(rho(:,end-6:end,4))*Mp)+B(:,end-6:end).*(-(rho(:,end-6:end,4))'*Mf)'-C(:,end-6:end).*rho(:,end-6:end,4);           
            Dur = A(:,end-6:end).*(-(u(:,end-6:end,4))*Mp)+B(:,end-6:end).*(-(u(:,end-6:end,4))'*Mf)'-C(:,end-6:end).*u(:,end-6:end,4);
            Dvr = A(:,end-6:end).*(-(v(:,end-6:end,4))*Mp)+B(:,end-6:end).*(-(v(:,end-6:end,4))'*Mf)'-C(:,end-6:end).*v(:,end-6:end,4);
            Dpr = A(:,end-6:end).*(-(p(:,end-6:end,4))*Mp)+B(:,end-6:end).*(-(p(:,end-6:end,4))'*Mf)'-C(:,end-6:end).*p(:,end-6:end,4);
        else %outflow boundary conditions
            Drhor = ((-Mx*rho(:,end-6:end,4)+Mx*p(:,end-6:end,4))*Mp)+A(:,end-6:end).*(-p(:,end-6:end,4)*Mp)+(((-My*rho(:,end-6:end,4)+My*p(:,end-6:end,4)))'*Mf)'+B(:,end-6:end).*(-p(:,end-6:end,4)'*Mf)'-C(:,end-6:end).*p(:,end-6:end,4);
            Dur = (-(Mx*u(:,end-6:end,4)+p(:,end-6:end,4))*Mp)+(-(My*u(:,end-6:end,4))'*Mf)';
            Dvr = (-(Mx*v(:,end-6:end,4))*Mp)+(-(My*v(:,end-6:end,4)+p(:,end-6:end,4))'*Mf)';
            Dpr = A(:,end-6:end).*(-(p(:,end-6:end,4))*Mp)+B(:,end-6:end).*(-(p(:,end-6:end,4))'*Mf)'-C(:,end-6:end).*p(:,end-6:end,4);
        end  
        
        %bottom boundary (radiation)
        Drhob = A(1:7,:).*(-(rho(1:7,:,4))*Mf)+B(1:7,:).*(-(rho(1:7,:,4))'*Mp)'-C(1:7,:).*rho(1:7,:,4);      
        Dub = A(1:7,:).*(-(u(1:7,:,4))*Mf)+B(1:7,:).*(-(u(1:7,:,4))'*Mp)'-C(1:7,:).*u(1:7,:,4);
        Dvb = A(1:7,:).*(-(v(1:7,:,4))*Mf)+B(1:7,:).*(-(v(1:7,:,4))'*Mp)'-C(1:7,:).*v(1:7,:,4);
        Dpb = A(1:7,:).*(-(p(1:7,:,4))*Mf)+B(1:7,:).*(-(p(1:7,:,4))'*Mp)'-C(1:7,:).*p(1:7,:,4);
        
        %top boundary
        if My == 0 %radiation boundary conditions
            Drhot = A(end-6:end,:).*(-(rho(end-6:end,:,4))*Mf)+B(end-6:end,:).*(-(rho(end-6:end,:,4))'*Mp)'-C(end-6:end,:).*rho(end-6:end,:,4);
            Dut = A(end-6:end,:).*(-(u(end-6:end,:,4))*Mf)+B(end-6:end,:).*(-(u(end-6:end,:,4))'*Mp)'-C(end-6:end,:).*u(end-6:end,:,4);
            Dvt = A(end-6:end,:).*(-(v(end-6:end,:,4))*Mf)+B(end-6:end,:).*(-(v(end-6:end,:,4))'*Mp)'-C(end-6:end,:).*v(end-6:end,:,4);
            Dpt = A(end-6:end,:).*(-(p(end-6:end,:,4))*Mf)+B(end-6:end,:).*(-(p(end-6:end,:,4))'*Mp)'-C(end-6:end,:).*p(end-6:end,:,4);
        else %outflow boundary conditions          
            Drhot = ((-Mx*rho(end-6:end,:,4)+Mx*p(end-6:end,:,4))*Mf)+A(end-6:end,:).*(-p(end-6:end,:,4)*Mf)+((-My*rho(end-6:end,:,4)+My*p(end-6:end,:,4))'*Mp)'+B(end-6:end,:).*(-p(end-6:end,:,4)'*Mp)'-C(end-6:end,:).*p(end-6:end,:,4);
            Dut = (-(Mx*u(end-6:end,:,4)+p(end-6:end,:,4))*Mf)+(-(My*u(end-6:end,:,4))'*Mp)';
            Dvt = (-(Mx*v(end-6:end,:,4))*Mf)+(-(My*v(end-6:end,:,4)+p(end-6:end,:,4))'*Mp)';
            Dpt = A(end-6:end,:).*(-(p(end-6:end,:,4))*Mf)+B(end-6:end,:).*(-(p(end-6:end,:,4))'*Mp)'-C(end-6:end,:).*p(end-6:end,:,4);
        end
        
        %stitch together interior and boundary nodes
        Drho(:,:,4) = [Drhol(:,1:3) [Drhob(1:3,4:end-3); Drhoi(4:end-3,4:end-3); Drhot(5:7,4:end-3)] Drhor(:,5:end)];
        Du(:,:,4) = [Dul(:,1:3) [Dub(1:3,4:end-3); Dui(4:end-3,4:end-3); Dut(5:7,4:end-3)] Dur(:,5:end)];
        Dv(:,:,4) = [Dvl(:,1:3) [Dvb(1:3,4:end-3); Dvi(4:end-3,4:end-3); Dvt(5:7,4:end-3)] Dvr(:,5:end)];
        Dp(:,:,4) = [Dpl(:,1:3) [Dpb(1:3,4:end-3); Dpi(4:end-3,4:end-3); Dpt(5:7,4:end-3)] Dpr(:,5:end)];        
        
        %calculate physical variables at new time period
        rhoswitch = rho(:,:,4)+dt*(b(1)*Drho(:,:,4)+b(2)*Drho(:,:,3)+b(3)*Drho(:,:,2)+b(4)*Drho(:,:,1));
        uswitch = u(:,:,4)+dt*(b(1)*Du(:,:,4)+b(2)*Du(:,:,3)+b(3)*Du(:,:,2)+b(4)*Du(:,:,1));
        vswitch = v(:,:,4)+dt*(b(1)*Dv(:,:,4)+b(2)*Dv(:,:,3)+b(3)*Dv(:,:,2)+b(4)*Dv(:,:,1));
        pswitch = p(:,:,4)+dt*(b(1)*Dp(:,:,4)+b(2)*Dp(:,:,3)+b(3)*Dp(:,:,2)+b(4)*Dp(:,:,1));
        
        %Save data
        Field(:,:,i) = pswitch(iField{1},iField{2});
        Sig(:,i) = pswitch(iSig{1},iSig{2});
   
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
            figure;pcolor(X,Y,rho(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Density contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            
            figure;pcolor(X,Y,u(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['U velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            
            figure;pcolor(X,Y,v(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['V velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
            
            figure;pcolor(X,Y,p(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Pressure contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
        end
    end
    save('out_pressure','Field','Sig','locations');
    
    compute_time = toc
end
