function [pomg,palphr,palphi, zreal, zimag, ysave] = comp2(m,nmode, omgall,alph,Rt,tr,omega)
%Inputs: Jet Mach number u_j/a_infity, nmode azimuthal mode, frequency range (omg), initial guessed alph 
%   This code is about inviscid compressible disturbance solution
%   search for complex eigenvalue alpha for specified range of real omg
%   For axisymmetric jet flow 
%==========summary of computation process==================================
%   There are a range of real omg(s)
%   For every omg, search for the right alpha using prediction-and-try method
%   For every guessed alpha, use RKF45 meathod to integrate from y=y1 to y=y2
%   For axisymmetric jet profiles
%==========summary of computation process==================================
%   Suggested Frequency Range (omg)
%   For mode 1: omgall=0.1:0.01:1.5
%   For mode 0: omgall=0.5:0.01:1.3
%   Suggested Initial guessed alph
%   For mode 1: alph=complex(0.13,-0.013)
%   For mode 0: alph=complex(0.6,-0.1)
% data ================================================
    
    num=length(omgall);
    % computation domain
    y1=0.005;
    y2=6;

    % error controls
    nmax=200;
    tol=1.e-5;
    abserr=1.e-5;
    % data for specific omg to be save for plot
    nomg=find(abs(omgall - omega) < 10*eps);%save eigenfunction data for omega = 0.85
    PYN=0;          % PYN=1 if data found for this omg
    % data ================================================
    % Loop of omg ############################################
    % Matrics to save true alph and corresponding omg
    talph=zeros(1,num);
    tomg=zeros(1,num);
    nn=0;  % number of omg whose right alph is found
    for kk=1:num
        omg=omgall(kk);
        %disp(['omg=',num2str(omg),', process=',num2str(kk),'/',num2str(num)])
        % extrapolate guessed starting alph for this omg
        if kk>2 %SETS THE INITIAL GUESS FOR ALPHA BASED ON THE PREVIOUSLY DETERMINED VALUES, USED PURELY FOR COMPUTATIONAL SPEED ?
            galphr=interp1(tomg(1:nn),real(talph(1:nn)),omg,'spline','extrap');
            galphi=interp1(tomg(1:nn),imag(talph(1:nn)),omg,'spline','extrap');
            alph=complex(galphr,galphi);
        end
        % Look for alph for this omg
        for itr=1:nmax+1
            % use guessed alph for c
            c=omg/alph;
            [u,u1]=u_sub(y1,Rt);
            [tep,tep1]=t_sub(m,u,u1,tr);
            umc=u-c;
            rkk=alph*alph*(1.-m*m*umc*umc/tep);
            wl=sqrt(rkk);
            dwla=((tep-m*m*umc*umc)*alph-m*m*omg*umc)/(wl*tep);
            % Impose Dirichlet BC at y=y1
            z1=zeros(4,1);
            z1(1)= besselj(nmode,wl*y1);
            z1(2)= wl*(besselj(nmode-1,wl*y1)-besselj(nmode+1,wl*y1))/2.;
            z1(3)= z1(2)*y1*dwla/wl;
            z1(4)= z1(3)/y1+wl*(besselj(nmode-2,wl*y1)-2.*besselj(nmode,wl*y1)+besselj(nmode+2,wl*y1))*y1*dwla/4.;
            % intergrate from y1 to y2, output (z,y) at every integration point
            [zz,yy]=RKF45(z1,y1,y2,abserr,alph,c,m,nmode,Rt,tr);
            % define variables
            [u,u1]=u_sub(y2,Rt);
            [tep,tep1]=t_sub(m,u,u1,tr);
            vmc=u-c;
            rll=alph*alph*(1.-m*m*vmc*vmc/tep);
            wm=sqrt(rll);
            % check accuracy of z at y=y2, guess new alph
            %er=besselh(nmode-1,wm*y2)-(nmode/wm/y2*zz(1,end)+zz(2,end)/wm);
            er=wm*zz(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-zz(2,end)*besselh(nmode,wm*y2);
            dwma=((tep-m*m*vmc*vmc)*alph-m*m*omg*vmc)/(wm*tep);
            %der=(besselh(nmode-2,y2*wm)-besselh(nmode,y2*wm))/2.*dwma*y2-(-nmode*zz(1,end)/(wm)^2/y2*dwma+nmode*zz(3,end)/wm/y2+zz(4,end)/wm-zz(2,end)/wm^2*dwma);
            der=dwma*zz(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.+wm*zz(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.+...
                wm*zz(1,end)*y2*dwma*(besselh(nmode-2,wm*y2)-2.*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4.-zz(4,end)*besselh(nmode,wm*y2)-zz(2,end)*dwma*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.;
            dalph=-er/der;

            alpherr=abs(dalph)/abs(alph);
            %disp(['itr=',' ',num2str(itr),', ','error=',num2str(alpherr),', abs error=',num2str(abs(zz(1,end)))])
            if (alpherr<tol)
                % right alph found, end search
                omgall(kk)
                break
            else
                % new alph for next trial
                alph=alph+dalph;
            end
        end
    

        % right alph found or not
        if itr>nmax
            % right alph not found, break loop of omg
            disp(['search failed: No alph found for this omg in ',num2str(nmax),' trials']);disp(' ')
            break
        else
            % right alph found
            nn=kk;
            %disp(['alph for current omg is: ', num2str(alph)])
            %disp(['# of trials to find alph for this omg= ', num2str(itr)]); disp(' ')
            talph(kk)=alph;
            tomg(kk)=omg;
            % For specific omg, save z (if found)
            if kk==nomg
                PYN=1;
                omgsave=omg;
                alphsave=alph;
                zsave=zz(1,:);
                ysave=yy;
            end
        end
    end
    % Loop of omg ############################################
    
    % plot alph against omg
    if nn>0
        % alph found for at least one omg
        % plot alphr against omg
        alphr=real(talph);
        alphi=imag(talph);
        figure(1)
        plot(tomg(1:nn),alphr(1:nn));
    %    axis([0,1,0,2])
        title('\alphar for a range of \omega')
        xlabel('\omega')
        ylabel('\alphar')
        grid on
    %   save alphr against omg
        pomg=tomg(1:nn);
        palphr=alphr(1:nn);
     %   save('omg_alphr','pomg','palphr')
     save omg_alphr pomg palphr;
     % plot alphi against omg
        figure(2)
        plot(tomg(1:nn),-alphi(1:nn));
     %   axis([0,1,0,0.4])
        title('-\alphai for a range of \omega')
        xlabel('\omega')
        ylabel('-\alphai')
        grid on    
    %   save alphr against omg
        palphi=-alphi(1:nn);
        save('omg_alphi','pomg','palphi')
    else
        disp('No alph found for any omg')
    end


    % plot eigenfunction for the specific omg
    if PYN==1
        % plot eigenfunction (complex form)
        figure(3)
        plot(zsave)
        title(['Eigenfunction (disturbance v) when \omg= ',num2str(omgsave)])
        xlabel('\phir') ; ylabel('\phii')
        grid on    
        % plot real part of eigenfunction against y
        figure(4)
        zreal=real(zsave);
        plot(zreal,ysave)
        title('real part of eigenfunction')
        xlabel('\phir'); ylabel('y');
        grid on    
        % plot imaginary part of eigenfunction against y
        figure(5)
        zimag=imag(zsave);
        plot(zimag,ysave)
        title('imaginary part of eigenfunction')
        xlabel('\phii'); ylabel('y');
        grid on    
        disp(' ')
        disp(['alph for omg [',num2str(omgsave),'] is [',num2str(alphsave),']']);disp(' ')
    else
        % alph for this specific omg not found
        disp(' ')
        omgspec=omgall(nomg);
        disp(['No alph found for the specific omg [',num2str(omgspec),'], No data to plot']);disp(' ')
    end

end

function [z,y]=RKF45(z1,ystart,yend,abserr,alpha,c,m,nmode,Rt,tr)
%   this code integrates a system of first order ordinary differential
%   equations by runge-kutta-fehlberg-45 method with automatic estimation
%   of local error and step size adjustment.
%   number of equations in the ODE system = size of z1.
%
%   [input]
%   z1 - vector of initial values of the dependent variables.
%   ystart - start value of the independent variable.
%   yend - end value of the independent variable.
%       When ystart>yend, use h= -0.05; while y(j)>yend; if y(j)+h<yend;
%       When ystart<yend, use h=  0.05; while y(j)<yend; if y(j)+h>yend;
%   abserr - bound of local error permitted, the computed solution
%       znew obtained in a step h must pass the test to be accepted
%       err<abserr
%   alpha - to be used in subfunction deriv, comes from the main program
%   c - to be used in subfunction deriv, comes from the main program
%
%   [output]
%   z - matrix of dependent variables at every integration point,
%       a column of z corresponds to a component in the y vector
%   y - vector of independent variable at every integration point.
%
%   [Comments]
%(1) The part of checking accuray and finding new h comes from document 10/16/2006
%(2) It's possible that: scheme for new h can not fulfill the error tolerance
%(3) To avoid dead loop, set maximun number of valid/total integrations.
%    When reaching maximun number of integration, go out of the function

    % parameters of RKF45 method
    a2=1/4; b2=1/4;
    a3=3/8; b3=3/32; c3=9/32;
    a4=12/13; b4=1932/2197; c4=-7200/2197; d4=7296/2197;
    a5=1; b5=439/216; c5=-8; d5=3680/513; e5=-845/4104;
    a6=1/2; b6=-8/27; c6=2; d6=-3544/2565; e6=1859/4104; f6=-11/40;
    n1=25/216; n3=1408/2565; n4=2197/4104; n5=-1/5;
    r1=1/360; r3=-128/4275; r4=-2197/75240; r5=1/50; r6=2/55;
    
    % step size
    h= 0.05;      % h can be either positive or negative
    hmin=0.001;   % hmin must be positive
    hmax=0.05;    % hmax must be positive
    
    % control parameters to avoid dead loop
    max1=1000;    % maximun number of valid integrations
    max2=2000;    % maximun number of all integrations
    i=0;          % number of all integration (including not valid integrations)
    j=1;          % number of valid integrations
    
    % results to be output
    y=[ ];
    y=[y,ystart];
    z=[ ];
    z=[z, z1];

    while y(j)<yend
        if y(j)+h>yend
            h=yend-y(j);
        end
        k1=h*deriv(y(j),z(:,j), alpha,c,m,nmode,Rt,tr);
        k2=h*deriv(y(j)+a2*h, z(:,j)+b2*k1, alpha,c,m,nmode,Rt,tr);
        k3=h*deriv(y(j)+a3*h, z(:,j)+b3*k1+c3*k2, alpha, c,m,nmode,Rt,tr);
        k4=h*deriv(y(j)+a4*h, z(:,j)+b4*k1+c4*k2+d4*k3, alpha, c,m,nmode,Rt,tr);
        k5=h*deriv(y(j)+a5*h, z(:,j)+b5*k1+c5*k2+d5*k3+e5*k4, alpha, c,m,nmode,Rt,tr);
        k6=h*deriv(y(j)+a6*h, z(:,j)+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5, alpha, c,m,nmode,Rt,tr);
        % new z at new y point
        znew=z(:,j)+n1*k1+n3*k3+n4*k4+n5*k5;
        % local error
        err=abs(r1*k1+r3*k3+r4*k4+r5*k5+r6*k6);
        err=max(err);
        % if err<abserr, accept new z; if not, go on to try new h
        if err<abserr
            y=[y,y(j)+h];
            z=[z,znew];
            j=j+1;
        end
        % new h (always s>0)
        s=(abserr*abs(h)/(2*err))^0.25;
        if s<0.1
            s=0.1;
        elseif s>4.0
            s=4.0;
        end
        h=s*h;
        if abs(h)>hmax
            h=hmax*h/abs(h);
        elseif abs(h)<hmin
            h=hmin*h/abs(h);
        end
        % reach maximum number of integration, break
        i=i+1;
        if j==max1 || i==max2
            disp(' ');
            disp('error in RKF45: too many integrations needed');
            return
        end
    end
end
    
    %=================================
function dz=deriv(y,z,alpha,c,m,nmode,Rt,tr)
% y - independent variable
% z - dependent variable
% alpha - complex wave number
% c = omega/alpha
% dz = dz/dy at y
    dz=zeros(4,1);
    [u,u1]=u_sub(y,Rt);
    [tep,tep1]=t_sub(m,u,u1,tr);
    cmc=u-c;
    dz(1)=z(2);                        
    dz(2)=(2*u1/cmc-tep1/tep-1./y)*z(2)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(1);  
    dz(3)=z(4);						   
    dz(4)=(2*u1/cmc-tep1/tep-1./y)*z(4)+2*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(3)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(3)...
        -z(1)*2*m^2*cmc*c*alpha/tep-2*u1*c*z(2)/(cmc*cmc*alpha);
end

function [u,u1]=u_sub(y,Rt)
% y - independent variable
% Rt - R/theta (jet radius / displacement thickness)
% u - velocity at y
% u2 - second derivative of u at y
    a=0.25*Rt;
    u=0.5*(1+tanh(a*(1./y-y)));
    u1=0.5*(1-(tanh(a*(1./y-y)))^2)*a*(-1./y^2-1);
end


function [tep,tep1]=t_sub(m,u,u1,tr)
% tr=T_infinity/T_jet
    g=(1.4-1)*m*m/2;
    tep=tr+(1-tr+g)*u-g*u^2;
    tep1=(1-tr+g)*u1-2*g*u*u1;
%b=1.5;
%tr=1.5;
%tep=1+(tr-1)*exp(-log(2)*(y/b).^2);
%tep1=2*(y/b*b)*(-log(2))*exp(-log(2)*(y/b).^2);
end