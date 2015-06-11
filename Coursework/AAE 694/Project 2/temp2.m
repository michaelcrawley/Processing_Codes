function [] = temp2(dx,dt,Mx,My)
    tic;
    dy = dx;
    X = -100:dx:100;
    Y = -100:dy:100;
    [x, y] = meshgrid(X,Y);
    time_itr_max = 0.05*2000; 
    foldername = [num2str(length(X)),' x ',num2str(length(Y)), ' Mx', strrep(num2str(Mx),'.','_'), ' My',strrep(num2str(My),'.','_')];
    
    
    %Boundary condition coefficients
%     theta = atan2(y,x);
    r = sqrt(x.^2+y.^2);
%     V = Mx*cos(theta)+My*sin(theta)+sqrt(1-(Mx*sin(theta)-My*cos(theta)).^2);
%     A = V.*cos(theta);
%     B = V.*sin(theta);
    V = Mx*(x./r)+My*(y./r)+sqrt(1-(Mx*(y./r)-My*(x./r)).^2);
    A = V.*(x./r);
    B = V.*(y./r);
    C = V./(2*r); 
    
    %3rd order optimized time discretization coefficients (b0, b1, b2, b3)
    b = [2.3025580883830 -2.4910075998482 1.5743409331815 -0.38589142217162];
    plot_times = 0.05*[1 100:100:1200 1600 2000 12000];

    %itialize variables
    rho = zeros(length(x),length(y),4);
    u = rho; v = rho; p = rho;
    
    rho(:,:,4) = 0.01*exp(-log(2)*(x.^2+y.^2)/9)+0.001*exp(-log(2)*((x-67).^2+y.^2)/25);
    u(:,:,4) = 0.0004*y.*exp(-log(2)*((x-67).^2+y.^2)/25);
    v(:,:,4) = -0.0004*(x-67).*exp(-log(2)*((x-67).^2+y.^2)/25);
    p(:,:,4) = 0.01*exp(-log(2)*(x.^2+y.^2)/9);
    
    g = figure;
    if exist(foldername,'file') ~= 7
        mkdir(pwd,foldername);
    end
    cd(foldername);
    h = waitbar(0,'Calculating...');
    imax = length(0:dt:time_itr_max)-1;
    for i = 1:imax %time step loop
        waitbar(i/imax,h,['i = ',num2str(i)]);
        %density step - interior nodes
        rhoswitchi = rho(:,:,4);
        for k=1:4
            rhoswitchi = rhoswitchi+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*rho(:,:,5-k)+u(:,:,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*rho(:,:,5-k)+v(:,:,5-k))','-DRP')');
        end

        %density step - left boundary (radiation)
        rhoswitchl = rho(:,1:7,4);
        for k=1:4
            rhoswitchl = rhoswitchl+(dt*b(k))*(A(:,1:7).*NumericalDerivative(1,6,dx,-(rho(:,1:7,5-k)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(rho(:,1:7,5-k))','-DRP')'-C(:,1:7).*rho(:,1:7,5-k));
        end
        
        if Mx == 0
            %density step - right boundary (radiation)
            rhoswitchr = rho(:,end-6:end,4);
            for k = 1:4
               rhoswitchr = rhoswitchr+(dt*b(k))*(A(:,end-6:end).*NumericalDerivative(1,6,dx,-(rho(:,end-6:end,5-k)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(rho(:,end-6:end,5-k))','-DRP')'-C(:,end-6:end).*rho(:,end-6:end,5-k));
            end            
        else
            %density step - right boundary (outflow)
            rhoswitchr = rho(:,end-6:end,4);
            for k = 1:4
               rhoswitchr = rhoswitchr+(dt*b(k))*(NumericalDerivative(1,6,dx,-Mx*rho(:,end-6:end,5-k)+Mx*p(:,end-6:end,5-k),'-DRP')+A(:,end-6:end).*NumericalDerivative(1,6,dx,-p(:,end-6:end,5-k),'-DRP')+NumericalDerivative(1,6,dy,(-My*rho(:,end-6:end,5-k)+My*p(:,end-6:end,5-k))','-DRP')'+B(:,end-6:end).*NumericalDerivative(1,6,dy,-p(:,end-6:end,5-k)','-DRP')'-C(:,end-6:end).*p(:,end-6:end,5-k));
            end
        end
        %density step - bottom boundary (radiation)
        rhoswitchb = rho(end-6:end,:,4);
        for k=1:4
            rhoswitchb = rhoswitchb +(dt*b(k))*(A(end-6:end,:).*NumericalDerivative(1,6,dx,-(rho(end-6:end,:,5-k)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(rho(end-6:end,:,5-k))','-DRP')'-C(end-6:end,:).*rho(end-6:end,:,5-k));
        end
        if My == 0
            %density step - top boundary (radiation)
            rhoswitcht = rho(1:7,:,4);
            for k = 1:4
               rhoswitcht = rhoswitcht+(dt*b(k))*(A(1:7,:).*NumericalDerivative(1,6,dx,-(rho(1:7,:,5-k)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(rho(1:7,:,5-k))','-DRP')'-C(1:7,:).*rho(1:7,:,5-k));
            end
        else
            %density step - top boundary (outflow)
            rhoswitcht = rho(1:7,:,4);
            for k = 1:4
               rhoswitcht = rhoswitcht+(dt*b(k))*(NumericalDerivative(1,6,dx,-Mx*rho(1:7,:,5-k)+Mx*p(1:7,:,5-k),'-DRP')+A(1:7,:).*NumericalDerivative(1,6,dx,-p(1:7,:,5-k),'-DRP')+NumericalDerivative(1,6,dy,(-My*rho(1:7,:,5-k)+My*p(1:7,:,5-k))','-DRP')'+B(1:7,:).*NumericalDerivative(1,6,dy,-p(1:7,:,5-k)','-DRP')'-C(1:7,:).*p(1:7,:,5-k));
            end
        end                                    
        
        rhonew = [rhoswitcht(1:3,:);rhoswitchl(4:end-3,1:3) rhoswitchi(4:end-3,4:end-3) rhoswitchr(4:end-3,5:end); rhoswitchb(5:end,:)];  

        
        %u velocity step - interior nodes
        uswitchi = u(:,:,4);
        for k =1:4
           uswitchi = uswitchi+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*u(:,:,5-k)+p(:,:,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*u(:,:,5-k))','-DRP')');
        end
        %u velocity step - left boundary (radiation)
        uswitchl = u(:,1:7,4);
        for k = 1:4
            uswitchl = uswitchl+(dt*b(k))*(A(:,1:7).*NumericalDerivative(1,6,dx,-(u(:,1:7,5-k)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(u(:,1:7,5-k))','-DRP')'-C(:,1:7).*u(:,1:7,5-k));
        end
        if Mx == 0
            %u velocity step - right boundary (radiation)
            uswitchr = u(:,end-6:end,4);
            for k = 1:4
               uswitchr = uswitchr+(dt*b(k))*(A(:,end-6:end).*NumericalDerivative(1,6,dx,-(u(:,end-6:end,5-k)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(u(:,end-6:end,5-k))','-DRP')'-C(:,end-6:end).*u(:,end-6:end,5-k));
            end 
        else
            %u velocity step - right boundary (outflow)
            uswitchr = u(:,end-6:end,4);
            for k = 1:4
               uswitchr = uswitchr+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*u(:,end-6:end,5-k)+p(:,end-6:end,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*u(:,end-6:end,5-k))','-DRP')');
            end
        end
        %u step - bottom boundary (radiation)
        uswitchb = u(end-6:end,:,4);
        for k=1:4
            uswitchb = uswitchb +(dt*b(k))*(A(end-6:end,:).*NumericalDerivative(1,6,dx,-(u(end-6:end,:,5-k)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(u(end-6:end,:,5-k))','-DRP')'-C(end-6:end,:).*u(end-6:end,:,5-k));
        end
        if My == 0
            %u velocity step - top boundary (radiation)
            uswitcht = u(1:7,:,4);
            for k = 1:4
               uswitcht = uswitcht+(dt*b(k))*(A(1:7,:).*NumericalDerivative(1,6,dx,-(u(1:7,:,5-k)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(u(1:7,:,5-k))','-DRP')'-C(1:7,:).*u(1:7,:,5-k));
            end
        else
            %u velocity step - top boundary (outflow)
            uswitcht = u(1:7,:,4);
            for k = 1:4
               uswitcht = uswitcht+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*u(1:7,:,5-k)+p(1:7,:,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*u(1:7,:,5-k))','-DRP')');
            end
        end
        unew = [uswitcht(1:3,:);uswitchl(4:end-3,1:3) uswitchi(4:end-3,4:end-3) uswitchr(4:end-3,5:end); uswitchb(5:end,:)];

        
        %v velocity step - interior nodes
        vswitchi = v(:,:,4);
        for k =1:4
           vswitchi = vswitchi+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*v(:,:,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*v(:,:,5-k)+p(:,:,5-k))','-DRP')');
        end
        %v velocity step - left boundary (radiation)
        vswitchl = v(:,1:7,4);
        for k = 1:4
            vswitchl = vswitchl+(dt*b(k))*(A(:,1:7).*NumericalDerivative(1,6,dx,-(v(:,1:7,5-k)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(v(:,1:7,5-k))','-DRP')'-C(:,1:7).*v(:,1:7,5-k));
        end
        if Mx == 0
            %v velocity step - right boundary (radiation)
            vswitchr = v(:,end-6:end,4);
            for k = 1:4
               vswitchr = vswitchr+(dt*b(k))*(A(:,end-6:end).*NumericalDerivative(1,6,dx,-(v(:,end-6:end,5-k)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(v(:,end-6:end,5-k))','-DRP')'-C(:,end-6:end).*v(:,end-6:end,5-k));
            end 
        else
            %v velocity step - right boundary (outflow)
            vswitchr = v(:,end-6:end,4);
            for k = 1:4
               vswitchr = vswitchr+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*v(:,end-6:end,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*v(:,end-6:end,5-k)+p(:,end-6:end,5-k))','-DRP')');
            end
        end
        %v step - bottom boundary (radiation)
        vswitchb = v(end-6:end,:,4);
        for k=1:4
            vswitchb = vswitchb +(dt*b(k))*(A(end-6:end,:).*NumericalDerivative(1,6,dx,-(v(end-6:end,:,5-k)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(v(end-6:end,:,5-k))','-DRP')'-C(end-6:end,:).*v(end-6:end,:,5-k));
        end
        if My == 0
            %v velocity step - top boundary (radiation)
            vswitcht = v(1:7,:,4);
            for k = 1:4
               vswitcht = vswitcht+(dt*b(k))*(A(1:7,:).*NumericalDerivative(1,6,dx,-(v(1:7,:,5-k)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(v(1:7,:,5-k))','-DRP')'-C(1:7,:).*v(1:7,:,5-k));
            end
        else
            %u velocity step - top boundary (outflow)
            vswitcht = v(1:7,:,4);
            for k = 1:4
               vswitcht = vswitcht+(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*v(1:7,:,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*v(1:7,:,5-k)+p(1:7,:,5-k))','-DRP')');
            end
        end        
        vnew = [vswitcht(1:3,:);vswitchl(4:end-3,1:3) vswitchi(4:end-3,4:end-3) vswitchr(4:end-3,5:end); vswitchb(5:end,:)];

        
        %pressure step - interior nodes
        pswitchi = p(:,:,4);
        for k = 1:4
           pswitchi = pswitchi +(dt*b(k))*(NumericalDerivative(1,6,dx,-(Mx*p(:,:,5-k)+u(:,:,5-k)),'-DRP')+NumericalDerivative(1,6,dy,-(My*p(:,:,5-k)+v(:,:,5-k))','-DRP')');
        end
        %pressure step - left boundary (radiation)
        pswitchl = p(:,1:7,4);
        for k = 1:4
           pswitchl = pswitchl+(dt*b(k))*(A(:,1:7).*NumericalDerivative(1,6,dx,-(p(:,1:7,5-k)),'-DRP')+B(:,1:7).*NumericalDerivative(1,6,dy,-(p(:,1:7,5-k))','-DRP')'-C(:,1:7).*p(:,1:7,5-k));
        end
        %pressure step - right boundary (outflow)
        pswitchr = p(:,end-6:end,4);
        for k = 1:4
           pswitchr = pswitchr+(dt*b(k))*(A(:,end-6:end).*NumericalDerivative(1,6,dx,-(p(:,end-6:end,5-k)),'-DRP')+B(:,end-6:end).*NumericalDerivative(1,6,dy,-(p(:,end-6:end,5-k))','-DRP')'-C(:,end-6:end).*p(:,end-6:end,5-k));
        end
        %pressure step - bottom boundary (radiation)
        pswitchb = p(end-6:end,:,4);
        for k=1:4
            pswitchb = pswitchb +(dt*b(k))*(A(end-6:end,:).*NumericalDerivative(1,6,dx,-(p(end-6:end,:,5-k)),'-DRP')+B(end-6:end,:).*NumericalDerivative(1,6,dy,-(p(end-6:end,:,5-k))','-DRP')'-C(end-6:end,:).*p(end-6:end,:,5-k));
        end
        %pressure step - top boundary (radiation)
        pswitcht = p(1:7,:,4);
        for k = 1:4
           pswitcht = pswitcht+(dt*b(k))*(A(1:7,:).*NumericalDerivative(1,6,dx,-(p(1:7,:,5-k)),'-DRP')+B(1:7,:).*NumericalDerivative(1,6,dy,-(p(1:7,:,5-k))','-DRP')'-C(1:7,:).*p(1:7,:,5-k));
        end
        pnew = [pswitcht(1:3,:);pswitchl(4:end-3,1:3) pswitchi(4:end-3,4:end-3) pswitchr(4:end-3,5:end); pswitchb(5:end,:)];                            
        
        %shift matrix to advance time step
        rho(:,:,1:3) = rho(:,:,2:4); 
        rho(:,:,4) = rhonew;        
        u(:,:,1:3) = u(:,:,2:4);
        u(:,:,4) = unew;
        v(:,:,1:3) = v(:,:,2:4);
        v(:,:,4) = vnew;       
        p(:,:,1:3) = p(:,:,2:4);
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
    cd ..
    close(h);
    ctime = toc
end