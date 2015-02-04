function [grab,X,Y,t] = Work4Cat5Prob1v2(dx)
%Linearized Euler Equation solver for NASA ICASE_LaRC CAA Workshop 4, benchmark Category 5, Problem 1 
%Version 2 uses dimensional variables, unlike version 1
%Inputs:    dx: non-dimensionalized grid spacing (axial and radial spacings will be equal)
    
    %Constants
    h = 0;
    b = 1.3;
    Rhalf = h+b;
    T_amb = 300;
    Tj = 600; %Remember to change this back to 360 after the code is validated
    gamma = 1.4;
    R = 287;
    a = sqrt(gamma*R*T_amb);
    Mj = 0.756;
    uj = Mj*sqrt(gamma*R*Tj);
    p_bar = 103330;
    rhoj = p_bar/R/Tj;
    rho_amb = p_bar/R/T_amb;
    M = uj/a;
    omega_o = 76; 
    Amp = 0.001; %Forcing amplitude
    Bx = 0.04*log(2);
    By = 0.32*log(2); 
    
    %Set up grid
    x = (-50:dx:150)*Rhalf;
    y = (0:dx:50)*Rhalf;
    [X,Y] = meshgrid(x,y);
    forcing = @(t) Amp*exp(-(Bx*X.^2 + By*Y.^2))*cos(omega_o*t); %forcing function
    
    %Find time-step
    bj = [2.3025580883830 -2.4910075998482 1.5743409331815 -0.38589142217162];
    Omega = timestepping(bj,M,length(x),0.1); %Fuck Tam
    dt = Omega*dx*Rhalf/(1.75*a*(M + sqrt(2))); 
    
    %Set maximum time for simulation
    imax = 175001;
    downsample = 500;
    dnsmax = (imax-1)/downsample+1;
    t = (0:imax-1)*dt;
    
    %Mean values    
    u_bar = (uj*exp(-log(2)*(Y/(b/Rhalf)-h/b).^2).*(Y > h) + uj.*(Y <= h));
    rho_bar = 1./(-0.5*(gamma-1)/(gamma*p_bar)*(u_bar-uj).*u_bar + 1/rhoj*u_bar/uj + 1/rho_amb * (uj-u_bar)/uj);
    
    %Normalized values
%     p_bar = p_bar/rhoj/a/a;
%     rho_bar = rho_bar/rhoj;
%     u_bar = u_bar/a;
    
    %Derivative matrices
    partial_x = mNumericalDerivative(1,6,dx,length(x),'-DRP'); %we'll have to transpose the order to get the multiplication right
    partial_y = mNumericalDerivative(1,6,dx,length(y),'-DRP');
    coefs = DRPlookup(0,6,dx)'; %asymmetric coefs for symmetric boundary condition
    dely_ubar = partial_y*u_bar;
    dely_rhobar = partial_y*rho_bar;

    %Boundary condition coefficients - NEED TO RECHECK THESE
    r = sqrt(X.^2+Y.^2);
    V = (u_bar/a).*(X./r)+sqrt(1-((u_bar/a).*(Y./r)).^2); %this should be Mach number, right?
    A = V.*(X./r);
    B = V.*(Y./r);
    C = V./(2*r); 
    
    %initialize variables - we'll leave these as initial conditions, since
    %the forcing should take care of everything.
    rho = zeros(length(y),length(x),1);
    u = rho; v = rho; 
    p = forcing(0); %start at t = 0
    Drho = zeros(length(y),length(x),4);
    Du = Drho;
    Dv = Drho;
    Dp = Drho;
    grab{1} = zeros(length(y),length(x),(imax-1)/downsample+1); %grab pressure perturbations 
    counter = 1;
    
    disp('Calculating....');
    for n = 1:imax        
        %time shift spatial derivative matrices
        Drho(:,:,1:3) = Drho(:,:,2:4);
        Du(:,:,1:3) = Du(:,:,2:4);
        Dv(:,:,1:3) = Dv(:,:,2:4);
        Dp(:,:,1:3) = Dp(:,:,2:4);
        
        %Partial derivatives
        delx_rho = (partial_x*rho')';
        dely_rho = partial_y*rho;
        delx_u = (partial_x*u')';
        dely_u = partial_y*u;        
        delx_v = (partial_x*v')';
        dely_v = partial_y*v;   
        delx_p = (partial_x*p')';
        dely_p = partial_y*p;
        
        %interior nodes
        Drho(:,:,4) = -u_bar.*delx_rho - rho_bar.*delx_u - dely_rhobar.*v - dely_v.*rho_bar;
        Du(:,:,4) = -u_bar.*delx_u - delx_p./rho_bar - dely_ubar.*v - dely_v.*u_bar;
        Dv(:,:,4) = -u_bar.*delx_v - dely_p./rho_bar;
        Dp(:,:,4) = forcing(t(n)) - gamma*p_bar.*delx_u - u_bar.*delx_p - gamma*p_bar*dely_v;
        
        %top boundary (radiation)        
        Drho(end-2:end,4:end-3,4) = -A(end-2:end,4:end-3).*delx_rho(end-2:end,4:end-3)-B(end-2:end,4:end-3).*dely_rho(end-2:end,4:end-3)-C(end-2:end,4:end-3).*rho(end-2:end,4:end-3);
        Du(end-2:end,4:end-3,4) = -A(end-2:end,4:end-3).*delx_u(end-2:end,4:end-3)-B(end-2:end,4:end-3).*dely_u(end-2:end,4:end-3)-C(end-2:end,4:end-3).*u(end-2:end,4:end-3);
        Dv(end-2:end,4:end-3,4) = -A(end-2:end,4:end-3).*delx_v(end-2:end,4:end-3)-B(end-2:end,4:end-3).*dely_v(end-2:end,4:end-3)-C(end-2:end,4:end-3).*v(end-2:end,4:end-3);
        Dp(end-2:end,4:end-3,4) = -A(end-2:end,4:end-3).*delx_p(end-2:end,4:end-3)-B(end-2:end,4:end-3).*dely_p(end-2:end,4:end-3)-C(end-2:end,4:end-3).*p(end-2:end,4:end-3);
        
        %right boundary (outflow)
        Dp(2:end,end-2:end,4) = -A(2:end,end-2:end).*delx_p(2:end,end-2:end)-B(2:end,end-2:end).*dely_p(2:end,end-2:end)-C(2:end,end-2:end).*p(2:end,end-2:end);
        Drho(:,end-2:end,4) = delx_p(:,end-2:end) - u_bar(:,end-2:end).*delx_rho(:,end-2:end) + dely_p(:,end-2:end)-A(:,end-2:end).*delx_p(:,end-2:end)-B(:,end-2:end).*dely_p(:,end-2:end)-C(:,end-2:end).*p(:,end-2:end);  
        
        %left boundary (radiation)
        Drho(:,1:3,4) = -A(:,1:3).*delx_rho(:,1:3)-B(:,1:3).*dely_rho(:,1:3)-C(:,1:3).*rho(:,1:3);
        Du(:,1:3,4) = -A(:,1:3).*delx_u(:,1:3)-B(:,1:3).*dely_u(:,1:3)-C(:,1:3).*u(:,1:3);
        Dv(:,1:3,4) = -A(:,1:3).*delx_v(:,1:3)-B(:,1:3).*dely_v(:,1:3)-C(:,1:3).*v(:,1:3);
        Dp(:,1:3,4) = -A(:,1:3).*delx_p(:,1:3)-B(:,1:3).*dely_p(:,1:3)-C(:,1:3).*p(:,1:3);  
        
        %Update variables
        rho = rho+dt*(bj(1)*Drho(:,:,4)+bj(2)*Drho(:,:,3)+bj(3)*Drho(:,:,2)+bj(4)*Drho(:,:,1));
        u = u+dt*(bj(1)*Du(:,:,4)+bj(2)*Du(:,:,3)+bj(3)*Du(:,:,2)+bj(4)*Du(:,:,1));
        v = v+dt*(bj(1)*Dv(:,:,4)+bj(2)*Dv(:,:,3)+bj(3)*Dv(:,:,2)+bj(4)*Dv(:,:,1));
        p = p+dt*(bj(1)*Dp(:,:,4)+bj(2)*Dp(:,:,3)+bj(3)*Dp(:,:,2)+bj(4)*Dp(:,:,1));
        
        %Set bottom boundary (symmetric)
        rho(1,:) = (-(1/coefs(1))*coefs(2:end)*rho(2:7,:));
        u(1,:) = (-(1/coefs(1))*coefs(2:end)*u(2:7,:));
        v(1,:) = 0;
        p(1,:) = (-(1/coefs(1))*coefs(2:end)*p(2:7,:));
        
        %Grab desired data
        if ~rem(n-1,downsample)
            disp([num2str(round(n/imax*10000)/100),'% Complete....']);
            grab{1}(:,:,counter) = p;
            counter = counter +1;
        end
    end
    t = t(1:downsample:end);
end

function [x0] = timestepping(bj,M,N,Delta)
%CAREFUL!!! YOU DON'T HAVE ANY ERROR CHECKING!!!

    %Equation 1 - indirect
    options = optimset('TolFun',1e-15,'MaxFunEvals',1000,'MaxIter',1000);
    dO = 0.0001;
    Omega_bar = (5*dO:dO:0.3)';
    soln = zeros(length(Omega_bar)+1,2);
    fval = zeros(length(Omega_bar),1);
    extflag = fval;

    soln(1,:) = [Omega_bar(1),-Omega_bar(1)/1000];
    for n = 1:length(Omega_bar)
        eqn = @(x) abs(1i*(exp(-1i*(x(1)+1i*x(2)))-1)/(bj(1) + bj(2)*exp(1i*(x(1)+1i*x(2))) + bj(3)*exp(1i*2*(x(1)+1i*x(2))) + bj(4)*exp(1i*3*(x(1)+1i*x(2)))) - Omega_bar(n));
        [soln(n+1,:),fval(n),extflag(n)] = fminsearch(eqn,soln(n,:),options);
    end

    soln = soln(2:end,:);
    Omega1 = soln(:,2);

    %Equation 2 - direct
    C = 1.75*(M + sqrt(2))*20*N/(Delta*log(10)*(1-M));
    Omega2 = -Omega_bar/C;

    %Find intersection
    [~,I] = min(abs(Omega2-Omega1));
    x0 = Omega_bar(I);
end