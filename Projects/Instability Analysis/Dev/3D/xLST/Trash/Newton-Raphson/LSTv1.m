function [talpha phi] = LSTv1(m,nmode,omega,alph,Rt,tr)
    %Inputs:
    %m: Mach number
    %nmode: Azimuthal mode
    %omgR: Omega range
    %alph: Initial guess for alpha
    %Rt: R/theta
    %tr: temperature ratio (T_infinity/T_jet)
    %version = T.1.0.0.0;
  

    %Velocity Profile
    U = @(y,Rt) 0.5*(1+tanh(0.25*Rt*(1./y-y)));
    %Temperature Profile
    T = @(y,Rt,tr) tr+(1-tr+((1.4-1)*m*m/2))*(U(y,Rt))-((1.4-1)*m*m/2)*(U(y,Rt))^2;
    
    y1 = .005;%Define range for velocity profile
    y2 = 4;    
    itrmax = 50; %Max iterations to find alpha
    tol=5.e-4; %tolerance for found alpha
    abserr=5.e-4;%tolerance for ode solver  
    
 

    for j = 1:itrmax
        alpha = alph;
        c = omega/alph;

        %Initial values for integration (along centerline of jet)
        u = U(y1,Rt);
        t = T(y1,Rt,tr);           
        umc=u-c;
        wl=sqrt(alph*alph*(1.-m*m*umc*umc/t));
        dwla=((t-m*m*umc*umc)*alph-m*m*omega*umc)/(wl*t);            
        z0(1,1)= besselj(nmode,wl*y1);
        z0(2,1)= wl*(besselj(nmode-1,wl*y1)-besselj(nmode+1,wl*y1))/2;
        z0(3,1)= z0(2,1)*y1*dwla/wl;
        z0(4,1)= z0(3,1)/y1+wl*(besselj(nmode-2,wl*y1)-2.*besselj(nmode,wl*y1)+besselj(nmode+2,wl*y1))*y1*dwla/4;

        %Boundary values for integration (at "inf")
        v = U(y2,Rt);
        t2 = T(y2,Rt,tr);            
        vmc = v-c;
        wm=sqrt(alph*alph*(1-m*m*vmc*vmc/t2));

        %Integration (centerline to inf)
        [z,y]=odeRKF45(@derive,z0,y1,y2,abserr,[alph omega m nmode Rt tr]); 

        %Calculate Errors            
        er=wm*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-z(2,end)*besselh(nmode,wm*y2);
        dwma=((t2-m*m*vmc*vmc)*alph-m*m*omega*vmc)/(wm*t2);
        der=dwma*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+wm*z(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+...
            wm*z(1,end)*y2*dwma*(besselh(nmode-2,wm*y2)-2*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4-z(4,end)*besselh(nmode,wm*y2)-z(2,end)*dwma*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2;
        dalpha=-er/der;
        alphaerr=abs(dalpha)/abs(alph);

        phi = z(1,:);
        psi = z(3,:);


        if (alphaerr<tol)% right alpha found, end search
            talpha = alpha;
            break            
        else% new alph for next trial
            alph=alph+dalpha;
            if imag(alph) > 0.1
                alph = complex(real(alph),-imag(alph));
            end
        end
        if j == itrmax %End program if convergence on alpha is not met
            disp('search failed: No alph found for this omg');
            talpha = NaN;
            pause(1);
            break
        end
    end
end

function dz=derive(y,z,inputs)
    % y - independent variable
    % z - dependent variable
    % alpha - complex wave number
    % c = omega/alpha
    % dz = dz/dy at y
    alpha = inputs(1);
    omega = inputs(2);
    m = inputs(3);
    nmode = inputs(4);
    Rt = inputs(5);
    tr = inputs(6);
    
    c = omega/alpha;
    
    %Velocity Profile
    a=0.25*Rt;
    u=0.5*(1+tanh(a*(1./y-y)));
    u1=0.5*(1-(tanh(a*(1./y-y)))^2)*a*(-1./y^2-1);
    %Temperature Profile
    g=((1.4-1)*m*m/2);
    tep=tr+(1-tr+g)*u-g*u^2;
    tep1=(1-tr+g)*u1-2*g*u*u1;
    
    dz=zeros(4,1);
    cmc=u-c;
    dz(1)=z(2);                        
    dz(2)=(2*u1/cmc-tep1/tep-1./y)*z(2)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(1);  
    dz(3)=z(4);						   
    dz(4)=(2*u1/cmc-tep1/tep-1./y)*z(4)+2*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(3)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(3)...
        -z(1)*2*m^2*cmc*c*alpha/tep-2*u1*c*z(2)/(cmc*cmc*alpha);
end

