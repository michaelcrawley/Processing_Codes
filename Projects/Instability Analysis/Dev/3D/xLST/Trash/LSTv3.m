function [alpha,phi] = LSTv3(m,nmode,tau,A,tr,omega,ialpha)
    %Analysis based on Spatial Linear Stability Theory, for Axisymmetric
    %Jet Flow.  This code uses curve fitted data.
    %inputs:
    %m: Mach Number
    %nmode: Azimuthal Mode
    %Rt: R/theta
    %tr: Temperature Ratio
    %omega: angular frequency
    %ialpha: initial guess for alpha

    y1 = .005;%Define range for velocity profile
    y2 = 4;    
    itrmax = 100; %Max iterations to find alpha
    tol=5e-4; %tolerance for found alpha
    abserr=5e-4;%tolerance for ode solver 
    
    U = @(y,tau,A)exp(-tau.*(y-A).^2).*(y>A)+1.*(y<=A);
    U1 =  @(y,tau,A)(tau.*(2*A - 2*y))./exp(tau.*(A - y).^2);
    T = @(y,tau,A,tr) tr+(1-tr+((1.4-1)*m*m/2))*(U(y,tau,A))-((1.4-1)*m*m/2)*(U(y,tau,A))^2;
    T1 = @(y,tau,A,m,tr) (1-tr+((1.4-1)*m*m/2))*U1(y,tau,A)-2*((1.4-1)*m*m/2)*U(y,tau,A)*U1(y,tau,A);    
    
    alphaerr = 1;    
    alpha = ialpha;
    counter = 1;
    while alphaerr >= tol && counter < itrmax
        c = omega/alpha;
        inputs = struct('alpha',alpha,'m',m,'nmode',nmode,'tau',tau,'A',A,'tr',tr,'c',c,'U',U,'U1',U1,'T',T,'T1',T1);
        %Initial values for integration (along centerline of jet)
        u = U(y1,tau,A);
        t = T(y1,tau,A,tr);           
        umc=u-c;
        wl=sqrt(alpha*alpha*(1.-m*m*umc*umc/t));
        dwla=((t-m*m*umc*umc)*alpha-m*m*omega*umc)/(wl*t);            
        z0(1,1)= besselj(nmode,wl*y1);
        z0(2,1)= wl*(besselj(nmode-1,wl*y1)-besselj(nmode+1,wl*y1))/2;
        z0(3,1)= z0(2,1)*y1*dwla/wl;
        z0(4,1)= z0(3,1)/y1+wl*(besselj(nmode-2,wl*y1)-2.*besselj(nmode,wl*y1)+besselj(nmode+2,wl*y1))*y1*dwla/4;
        
        %Boundary values for integration (at "inf")
        v = U(y2,tau,A);
        t2 = T(y2,tau,A,tr);            
        vmc = v-c;
        wm=sqrt(alpha*alpha*(1-m*m*vmc*vmc/t2));
        
        %Integration (centerline to inf)
        [z,y]=odeRKF45(@derive,z0,y1,y2,abserr,inputs); 
        phi = z(1,:);
        
        %Calculate Errors            
        er=wm*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-z(2,end)*besselh(nmode,wm*y2);
        dwma=((t2-m*m*vmc*vmc)*alpha-m*m*omega*vmc)/(wm*t2);
        der=dwma*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+wm*z(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+...
            wm*z(1,end)*y2*dwma*(besselh(nmode-2,wm*y2)-2*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4-z(4,end)*besselh(nmode,wm*y2)-z(2,end)*dwma*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2;
        dalpha=-er/der;
        alphaerr=abs(dalpha)/abs(alpha);
        
        if alphaerr>tol
            alpha=alpha+dalpha;
        elseif imag(alpha) > 0.1
            alpha = complex(real(alpha),-imag(alpha));
        end
        counter = counter +1;
    end
end

function dz=derive(y,z,inputs)
    % y - independent variable
    % z - dependent variable
    % alpha - complex wave number
    % c = omega/alpha
    % dz = dz/dy at y
    alpha = inputs.alpha;
    m = inputs.m;
    nmode = inputs.nmode;
    tau = inputs.tau;
    A = inputs.A;
    tr = inputs.tr;
    c = inputs.c;
    
    %Velocity Profile
    u=inputs.U(y,tau,A);
    u1=inputs.U1(y,tau,A);
    %Temperature Profile
    tep=inputs.T(y,tau,A,tr);
    tep1=inputs.T1(y,tau,A,m,tr);
    
    dz=zeros(4,1);
    cmc=u-c;
    dz(1)=z(2);                        
    dz(2)=(2*u1/cmc-tep1/tep-1./y)*z(2)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(1);  
    dz(3)=z(4);						   
    dz(4)=(2*u1/cmc-tep1/tep-1./y)*z(4)+2*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(3)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(3)...
        -z(1)*2*m^2*cmc*c*alpha/tep-2*u1*c*z(2)/(cmc*cmc*alpha);
end