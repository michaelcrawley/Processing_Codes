function [] = temp3(dx,dt,Mx,My)
    tic;
    dy = dx;
    X = -100:dx:100;
    Y = -100:dy:100;
    [x, y] = meshgrid(X,Y);
    time_itr_max = 0.05*2000; 
    foldername = [num2str(length(X)),' x ',num2str(length(Y)), ' Mx', strrep(num2str(Mx),'.','_'), ' My',strrep(num2str(My),'.','_')];
    
    
    %Boundary condition coefficients
    r = sqrt(x.^2+y.^2);
    V = Mx*(x./r)+My*(y./r)+sqrt(1-(Mx*(y./r)-My*(x./r)).^2);
    A = V.*(x./r);
    B = V.*(y./r);
    C = V./(2*r); 
    
    %3rd order optimized time discretization coefficients (b0, b1, b2, b3)
    b = [2.3025580883830 -2.4910075998482 1.5743409331815 -0.38589142217162];
    plot_times = 0.05*[1 100:100:1200 1600 2000 12000];

    %itialize variables
    rho = zeros(length(y),length(x),4);
    u = rho; v = rho; p = rho;

    shift = [0 0 -1];
    %rho spatial derivatives
    Drhoi = rho;
    Drhol = zeros(length(y),7,4);
    Drhor = Drhol;
    Drhob = zeros(7,length(x),4);
    Drhot = Drhob;
    
    %u velocity spatial derivatives
    Dui = u;
    Dul = Drhol;
    Dur = Drhol;
    Dub = Drhob;
    Dut = Drhob;
    
    %u velocity spatial derivatives
    Dvi = v;
    Dvl = Drhol;
    Dvr = Drhol;
    Dvb = Drhob;
    Dvt = Drhob;
    
    %pressure spatial derivatives
    Dpi = p;
    Dpl = Drhol;
    Dpr = Drhol;
    Dpb = Drhob;
    Dpt = Drhob;
    
    rho(:,:,4) = 0.01*exp(-log(2)*(x.^2+y.^2)/9)+0.001*exp(-log(2)*((x-67).^2+y.^2)/25);
    u(:,:,4) = 0.0004*y.*exp(-log(2)*((x-67).^2+y.^2)/25);
    v(:,:,4) = -0.0004*(x-67).*exp(-log(2)*((x-67).^2+y.^2)/25);
    p(:,:,4) = 0.01*exp(-log(2)*(x.^2+y.^2)/9);
    
    g = figure;
%     if exist(foldername,'file') ~= 7
%         mkdir(pwd,foldername);
%     end
%     cd(foldername);
    h = waitbar(0,'Calculating...');
    imax = length(0:dt:time_itr_max)-1;
    for i = 1:imax %time step loop
        waitbar(i/imax,h,['i = ',num2str(i)]);
        
        %time shift spatial derivative matrices
        Drhoi = circshift(Drhoi,shift);
        Drhol = circshift(Drhol,shift);
        Drhor = circshift(Drhor,shift);
        Drhob = circshift(Drhob,shift);
        Drhot = circshift(Drhot,shift);
        Dui = circshift(Dui,shift);
        Dul = circshift(Dul,shift);
        Dur = circshift(Dur,shift);
        Dub = circshift(Dub,shift);
        Dut = circshift(Dut,shift);
        Dvi = circshift(Dvi,shift);
        Dvl = circshift(Dvl,shift);
        Dvr = circshift(Dvr,shift);
        Dvb = circshift(Dvb,shift);
        Dvt = circshift(Dvt,shift);
        Dpi = circshift(Dpi,shift);
        Dpl = circshift(Dpl,shift);
        Dpr = circshift(Dpr,shift);
        Dpb = circshift(Dpb,shift);
        Dpt = circshift(Dpt,shift);
        
        %%Density
        %interior nodes
        Drhoi(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*rho(:,:,4)+u(:,:,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*rho(:,:,4)+v(:,:,4))','-DRP')';        
        rhoswitchi = rho(:,:,4)+dt*(b(1)*Drhoi(:,:,4)+b(2)*Drhoi(:,:,3)+b(3)*Drhoi(:,:,2)+b(4)*Drhoi(:,:,1));
        %left boundary (radiation)
        Drhol(:,:,4) = A(:,1:7).*NumericalDerivative(1,6,dx,-(rho(:,1:7,4)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(rho(:,1:7,4))','-DRP')'-C(:,1:7).*rho(:,1:7,4);
        rhoswitchl = rho(:,1:7,4)+dt*(b(1)*Drhol(:,:,4)+b(2)*Drhol(:,:,3)+b(3)*Drhol(:,:,2)+b(4)*Drhol(:,:,1));        
        %right boundary
        if Mx == 0 %radiation boundary conditions            
            Drhor(:,:,4) = A(:,end-6:end).*NumericalDerivative(1,6,dx,-(rho(:,end-6:end,4)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(rho(:,end-6:end,4))','-DRP')'-C(:,end-6:end).*rho(:,end-6:end,4);           
        else %outflow boundary conditions
            Drhor(:,:,4) = NumericalDerivative(1,6,dx,-Mx*rho(:,end-6:end,4)+Mx*p(:,end-6:end,4),'-DRP')+A(:,end-6:end).*NumericalDerivative(1,6,dx,-p(:,end-6:end,4),'-DRP')+NumericalDerivative(1,6,dy,(-My*rho(:,end-6:end,4)+My*p(:,end-6:end,4))','-DRP')'+B(:,end-6:end).*NumericalDerivative(1,6,dy,-p(:,end-6:end,4)','-DRP')'-C(:,end-6:end).*p(:,end-6:end,4);
        end
        rhoswitchr = rho(:,end-6:end,4)+dt*(b(1)*Drhor(:,:,4)+b(2)*Drhor(:,:,3)+b(3)*Drhor(:,:,2)+b(4)*Drhor(:,:,1));        
        %bottom boundary (radiation)
        Drhob(:,:,4) = A(end-6:end,:).*NumericalDerivative(1,6,dx,-(rho(end-6:end,:,4)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(rho(end-6:end,:,4))','-DRP')'-C(end-6:end,:).*rho(end-6:end,:,4);
        rhoswitchb = rho(end-6:end,:,4)+dt*(b(1)*Drhob(:,:,4)+b(2)*Drhob(:,:,3)+b(3)*Drhob(:,:,2)+b(4)*Drhob(:,:,1));        
        %top boundary
        if My == 0 %radiation boundary conditions
            Drhot(:,:,4) = A(1:7,:).*NumericalDerivative(1,6,dx,-(rho(1:7,:,4)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(rho(1:7,:,4))','-DRP')'-C(1:7,:).*rho(1:7,:,4);
        else %outflow boundary conditions          
            Drhot(:,:,4) = NumericalDerivative(1,6,dx,-Mx*rho(1:7,:,4)+Mx*p(1:7,:,4),'-DRP')+A(1:7,:).*NumericalDerivative(1,6,dx,-p(1:7,:,4),'-DRP')+NumericalDerivative(1,6,dy,(-My*rho(1:7,:,4)+My*p(1:7,:,4))','-DRP')'+B(1:7,:).*NumericalDerivative(1,6,dy,-p(1:7,:,4)','-DRP')'-C(1:7,:).*p(1:7,:,4);
        end
        rhoswitcht = rho(1:7,:,4)+dt*(b(1)*Drhot(:,:,4)+b(2)*Drhot(:,:,3)+b(3)*Drhot(:,:,2)+b(4)*Drhot(:,:,1));

        %%u velocity
        %interior nodes
        Dui(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*u(:,:,4)+p(:,:,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*u(:,:,4))','-DRP')';
        uswitchi = u(:,:,4)+dt*(b(1)*Dui(:,:,4)+b(2)*Dui(:,:,3)+b(3)*Dui(:,:,2)+b(4)*Dui(:,:,1));
        %left boundary (radiation)
        Dul(:,:,4) = A(:,1:7).*NumericalDerivative(1,6,dx,-(u(:,1:7,4)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(u(:,1:7,4))','-DRP')'-C(:,1:7).*u(:,1:7,4);
        uswitchl = u(:,1:7,4)+dt*(b(1)*Dul(:,:,4)+b(2)*Dul(:,:,3)+b(3)*Dul(:,:,2)+b(4)*Dul(:,:,1));        
        %right boundary
        if Mx == 0 %radiation boundary conditions            
            Dur(:,:,4) = A(:,end-6:end).*NumericalDerivative(1,6,dx,-(u(:,end-6:end,4)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(u(:,end-6:end,4))','-DRP')'-C(:,end-6:end).*u(:,end-6:end,4); 
        else %outflow boundary conditions
            Dur(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*u(:,end-6:end,4)+p(:,end-6:end,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*u(:,end-6:end,4))','-DRP')';
        end
        uswitchr = u(:,end-6:end,4)+dt*(b(1)*Dur(:,:,4)+b(2)*Dur(:,:,3)+b(3)*Dur(:,:,2)+b(4)*Dur(:,:,1));        
        %bottom boundary (radiation)
        Dub(:,:,4) = A(end-6:end,:).*NumericalDerivative(1,6,dx,-(u(end-6:end,:,4)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(u(end-6:end,:,4))','-DRP')'-C(end-6:end,:).*u(end-6:end,:,4);
        uswitchb = u(end-6:end,:,4)+dt*(b(1)*Dub(:,:,4)+b(2)*Dub(:,:,3)+b(3)*Dub(:,:,2)+b(4)*Dub(:,:,1));
        %top boundary
        if My == 0 %radiation boundary conditions
            Dut(:,:,4) = A(1:7,:).*NumericalDerivative(1,6,dx,-(u(1:7,:,4)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(u(1:7,:,4))','-DRP')'-C(1:7,:).*u(1:7,:,4);
        else %outflow boundary conditions
            Dut(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*u(1:7,:,4)+p(1:7,:,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*u(1:7,:,4))','-DRP')';
        end
        uswitcht = u(1:7,:,4)+dt*(b(1)*Dut(:,:,4)+b(2)*Dut(:,:,3)+b(3)*Dut(:,:,2)+b(4)*Dut(:,:,1));

        %%v velocity
        %interior nodes
        Dvi(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*v(:,:,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*v(:,:,4)+p(:,:,4))','-DRP')';
        vswitchi = v(:,:,4)+dt*(b(1)*Dvi(:,:,4)+b(2)*Dvi(:,:,3)+b(3)*Dvi(:,:,2)+b(4)*Dvi(:,:,1));        
        %left boundary (radiation)
        Dvl(:,:,4) = A(:,1:7).*NumericalDerivative(1,6,dx,-(v(:,1:7,4)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(v(:,1:7,4))','-DRP')'-C(:,1:7).*v(:,1:7,4);
        vswitchl = v(:,1:7,4)+dt*(b(1)*Dvl(:,:,4)+b(2)*Dvl(:,:,3)+b(3)*Dvl(:,:,2)+b(4)*Dvl(:,:,1));
        %right boundary
        if Mx == 0 %radiation boundary conditions
            Dvr(:,:,4) = A(:,end-6:end).*NumericalDerivative(1,6,dx,-(v(:,end-6:end,4)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(v(:,end-6:end,4))','-DRP')'-C(:,end-6:end).*v(:,end-6:end,4);
        else %outflow boundary conditions
            Dvr(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*v(:,end-6:end,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*v(:,end-6:end,4)+p(:,end-6:end,4))','-DRP')';
        end
        vswitchr = v(:,end-6:end,4)+dt*(b(1)*Dvr(:,:,4)+b(2)*Dvr(:,:,3)+b(3)*Dvr(:,:,2)+b(4)*Dvr(:,:,1));       
        %bottom boundary (radiation)
        Dvb(:,:,4) = A(end-6:end,:).*NumericalDerivative(1,6,dx,-(v(end-6:end,:,4)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(v(end-6:end,:,4))','-DRP')'-C(end-6:end,:).*v(end-6:end,:,4);
        vswitchb = v(end-6:end,:,4)+dt*(b(1)*Dvb(:,:,4)+b(2)*Dvb(:,:,3)+b(3)*Dvb(:,:,2)+b(4)*Dvb(:,:,1));
        %top boundary
        if My == 0 %radiation boundary conditions
            Dvt(:,:,4) = A(1:7,:).*NumericalDerivative(1,6,dx,-(v(1:7,:,4)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(v(1:7,:,4))','-DRP')'-C(1:7,:).*v(1:7,:,4);
        else %outflow boundary conditions
            Dvt(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*v(1:7,:,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*v(1:7,:,4)+p(1:7,:,4))','-DRP')';
        end    
        vswitcht = v(1:7,:,4)+dt*(b(1)*Dvt(:,:,4)+b(2)*Dvt(:,:,3)+b(3)*Dvt(:,:,2)+b(4)*Dvt(:,:,1));

        %%pressure
        %interior nodes
        Dpi(:,:,4) = NumericalDerivative(1,6,dx,-(Mx*p(:,:,4)+u(:,:,4)),'-DRP')+NumericalDerivative(1,6,dy,-(My*p(:,:,4)+v(:,:,4))','-DRP')';
        pswitchi = p(:,:,4)+dt*(b(1)*Dpi(:,:,4)+b(2)*Dpi(:,:,3)+b(3)*Dpi(:,:,2)+b(4)*Dpi(:,:,1));        
        %left boundary (radiation)
        Dpl(:,:,4) = A(:,1:7).*NumericalDerivative(1,6,dx,-(p(:,1:7,4)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(p(:,1:7,4))','-DRP')'-C(:,1:7).*p(:,1:7,4);
        pswitchl = p(:,1:7,4)+dt*(b(1)*Dpl(:,:,4)+b(2)*Dpl(:,:,3)+b(3)*Dpl(:,:,2)+b(4)*Dpl(:,:,1));        
        %right boundary (outflow)
        Dpr(:,:,4) = A(:,end-6:end).*NumericalDerivative(1,6,dx,-(p(:,end-6:end,4)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(p(:,end-6:end,4))','-DRP')'-C(:,end-6:end).*p(:,end-6:end,4);
        pswitchr = p(:,end-6:end,4)+dt*(b(1)*Dpr(:,:,4)+b(2)*Dpr(:,:,3)+b(3)*Dpr(:,:,2)+b(4)*Dpr(:,:,1));        
        %bottom boundary (radiation)
        Dpb(:,:,4) = A(end-6:end,:).*NumericalDerivative(1,6,dx,-(p(end-6:end,:,4)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(p(end-6:end,:,4))','-DRP')'-C(end-6:end,:).*p(end-6:end,:,4);
        pswitchb = p(end-6:end,:,4)+dt*(b(1)*Dpb(:,:,4)+b(2)*Dpb(:,:,3)+b(3)*Dpb(:,:,2)+b(4)*Dpb(:,:,1));        
        %top boundary (radiation)
        Dpt(:,:,4) = A(1:7,:).*NumericalDerivative(1,6,dx,-(p(1:7,:,4)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(p(1:7,:,4))','-DRP')'-C(1:7,:).*p(1:7,:,4);
        pswitcht = p(1:7,:,4)+dt*(b(1)*Dpt(:,:,4)+b(2)*Dpt(:,:,3)+b(3)*Dpt(:,:,2)+b(4)*Dpt(:,:,1));                                 

        %stitch boundary nodes and interior nodes together
        rhonew = [rhoswitchl(:,1:3) [rhoswitcht(1:3,4:end-3); rhoswitchi(4:end-3,4:end-3); rhoswitchb(5:7,4:end-3)] rhoswitchr(:,5:7)];
        unew = [uswitchl(:,1:3) [uswitcht(1:3,4:end-3); uswitchi(4:end-3,4:end-3); uswitchb(5:7,4:end-3)] uswitchr(:,5:7)];
        vnew = [vswitchl(:,1:3) [vswitcht(1:3,4:end-3); vswitchi(4:end-3,4:end-3); vswitchb(5:7,4:end-3)] vswitchr(:,5:7)];
        pnew = [pswitchl(:,1:3) [pswitcht(1:3,4:end-3); pswitchi(4:end-3,4:end-3); pswitchb(5:7,4:end-3)] pswitchr(:,5:7)];  
        
        
        %shift physical variable matrices
        rho = circshift(rho,shift); 
        u = circshift(u,shift);
        v = circshift(v,shift);
        p = circshift(p,shift);
        rho(:,:,4) = rhonew;        
        u(:,:,4) = unew;        
        v(:,:,4) = vnew;        
        p(:,:,4) = pnew;        
        
        figure(g);
        surf(X,Y,rho(:,:,4)); shading interp; view([45 45]);drawnow;title(['Density profile, i = ',num2str(i)]);
%         saveas(g,['densitysurf_i' num2str(i)]);

%         if any(abs(dt*i - plot_times) <= eps)
%             pcolor(X,Y,rho(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Density contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
%             saveas(gcf,['density_i',num2str(i)],'fig');
%             
%             pcolor(X,Y,u(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['U velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
%             saveas(gcf,['u_i',num2str(i)],'fig');
%             
%             pcolor(X,Y,v(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['V velocity contour for t = ',num2str(dt*i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
%             saveas(gcf,['v_i',num2str(i)],'fig');
%             
%             pcolor(X,Y,p(:,:,end)); colorbar; shading interp; colormap(flipud(gray)); xlabel('x'); ylabel('y'); title(['Pressure contour for t = ',num2str(i),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
%             saveas(gcf,['p_i',num2str(i)],'fig');
%             close all;
%         end
    end
%     cd ..
    close(h); close(g);
    ctime = toc
end