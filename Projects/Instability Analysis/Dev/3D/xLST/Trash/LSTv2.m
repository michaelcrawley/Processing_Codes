function [alpha,phi] = LSTv2(m,nmode,tr,Ufun,Uparams,U1fun,U1params,yPIV,omega,ialpha)
    %Analysis based on Spatial Linear Stability Theory, for Axisymmetric
    %Jet Flow. This code uses PIV data for the velocity profile.
    %inputs:
    %m: Mach Number
    %nmode: Azimuthal Mode
    %tr: Temperature Ratio
    %U: velocity profile (unscaled)
    %yPIV: radial positions for U data points (scaled)
    %omega: angular frequency
    %ialpha: initial guess for alpha

    y1 = .005;%Define range for velocity profile
    y2 = 4;    
    itrmax = 100; %Max iterations to find alpha
    tol=5e-4; %tolerance for found alpha
    abserr=5e-4;%tolerance for ode solver
    
    T = @(y,U,tr) tr+(1-tr+((1.4-1)*m*m/2))*(U)-((1.4-1)*m*m/2)*(U)^2;
    T1 = @(y,U,U1,tr) (1-tr+((1.4-1)*m*m/2))*U1-2*((1.4-1)*m*m/2)*U*U1;     
    
    u = Ufun(y1,Uparams);
    v = Ufun(y2,Uparams); 
    t = T(y1,u,tr);
    t2 = T(y2,v,tr);
    
    alphaerr = 1;    
    alpha = ialpha;
    counter = 1;
    while alphaerr >= tol && counter < itrmax
        c = omega/alpha;
        inputs = struct('alpha',alpha,'m',m,'nmode',nmode,'tr',tr,'c',c,'T',T,'T1',T1,'yPIV',yPIV,'Ufun',Ufun,'Uparams',Uparams,'U1fun',U1fun,'U1params',U1params);
        
        %Initial values for integration (along centerline of jet)          
        umc=u-c;
        wl=sqrt(alpha*alpha*(1.-m*m*umc*umc/t));
        dwla=((t-m*m*umc*umc)*alpha-m*m*omega*umc)/(wl*t);            
        z0(1,1)= besselj(nmode,wl*y1);
        z0(2,1)= wl*(besselj(nmode-1,wl*y1)-besselj(nmode+1,wl*y1))/2;
        z0(3,1)= z0(2,1)*y1*dwla/wl;
        z0(4,1)= z0(3,1)/y1+wl*(besselj(nmode-2,wl*y1)-2.*besselj(nmode,wl*y1)+besselj(nmode+2,wl*y1))*y1*dwla/4;
        
        %Boundary values for integration (at "inf")           
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
        if counter >= itrmax
           disp('No solution found');
           alpha = NaN;
        end
    end
end

function [dz]=derive(y,z,inputs)
    % y - independent variable
    % z - dependent variable
    % alpha - complex wave number
    % c = omega/alpha
    % dz = dz/dy at y
    alpha = inputs.alpha;
    m = inputs.m;
    nmode = inputs.nmode;
    c = inputs.c;
    tr = inputs.tr;
    
    %Velocity Profile
    u=inputs.Ufun(y,inputs.Uparams);
    u1=inputs.U1fun(y,inputs.U1params);
    %Temperature Profile
    tep=inputs.T(y,u,tr);
    tep1=inputs.T1(y,u,u1,tr);
    
    dz=zeros(4,1);
    cmc=u-c;
    dz(1)=z(2);                        
    dz(2)=(2*u1/cmc-tep1/tep-1./y)*z(2)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(1);  
    dz(3)=z(4);						   
    dz(4)=(2*u1/cmc-tep1/tep-1./y)*z(4)+2*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(3)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(3)...
        -z(1)*2*m^2*cmc*c*alpha/tep-2*u1*c*z(2)/(cmc*cmc*alpha);
end