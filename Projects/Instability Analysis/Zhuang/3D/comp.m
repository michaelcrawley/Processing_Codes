%   This code is about inviscid compressible disturbance solution
%   search for complex eigenvalue alpha for specified range of real omg
%   For axisymmetric jet flow 
%==========summary of computation process==================================
%   There are a range of real omg(s)
%   For every omg, search for the right alpha using prediction-and-try method
%   For every guessed alpha, use RKF45 meathod to integrate from y=y1 to y=y2
%   For axisymmetric jet profiles
%==========summary of computation process==================================

clear all
close all
pause(1)

% data ================================================
% frequency range (omg)
%omgall=0.1:0.01:1.5; %mode=1
omgall=0.5:0.01:1.3; %mode=0
num=length(omgall);
% initial guessed alph
%alph=complex(0.13,-0.013); %mode=1
alph=complex(0.6,-0.1); %mode=0
% computation domain
y1=0.005;
y2=6;
%Jet Mach number u_j/a_infity
m=0;
%nmode azimuthal mode
nmode=0;
% error controls
nmax=200;
tol=1.e-5;
abserr=1.e-5;
% data for specific omg to be save for plot
nomg=6;
PYN=0;          % PYN=1 if data found for this omg
% data ================================================
% Loop of omg ############################################
% Matrics to save true alph and coresponding omg
talph=zeros(1,num);
tomg=zeros(1,num);
nn=0;  % number of omg whose right alph is found
for kk=1:num
    omg=omgall(kk);
    disp(['omg=',num2str(omg),', process=',num2str(kk),'/',num2str(num)])
    % extrapolate guessed starting alph for this omg
    if kk>2
        galphr=interp1(tomg(1:nn),real(talph(1:nn)),omg,'spline','extrap');
        galphi=interp1(tomg(1:nn),imag(talph(1:nn)),omg,'spline','extrap');
        alph=complex(galphr,galphi);
    end
    % Look for alph for this omg
    for itr=1:nmax+1
        % use guessed alph for c
        c=omg/alph;
        [u,u1]=u_sub(y1);
        [tep,tep1]=t_sub(y1,m);
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
        [zz,yy]=RKF45(z1,y1,y2,abserr,alph,c,m,nmode);
        % define variables
        [u,u1]=u_sub(y2);
        [tep,tep1]=t_sub(y2,m);
        vmc=u-c;
        rll=alph*alph*(1.-m*m*vmc*vmc/tep);
        wm=sqrt(rll);
        % check accuracy of z at y=y2, guess new alph
%        er=besselh(nmode-1,wm*y2)-(nmode/wm/y2*zz(1,end)+zz(2,end)/wm);
        er=wm*zz(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-zz(2,end)*besselh(nmode,wm*y2)
        dwma=((tep-m*m*vmc*vmc)*alph-m*m*omg*vmc)/(wm*tep);
%        der=(besselh(nmode-2,y2*wm)-besselh(nmode,y2*wm))/2.*dwma*y2-(-nmode*zz(1,end)/(wm)^2/y2*dwma+nmode*zz(3,end)/wm/y2+zz(4,end)/wm-zz(2,end)/wm^2*dwma);
        der=dwma*zz(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.+wm*zz(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.+...
            wm*zz(1,end)*y2*dwma*(besselh(nmode-2,wm*y2)-2.*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4.-zz(4,end)*besselh(nmode,wm*y2)-zz(2,end)*dwma*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.
        dalph=-er/der;
%        alpherr=abs(dalph)/abs(alph);
        alpherr=abs(dalph)/abs(alph);
        disp(['itr=',' ',num2str(itr),', ','error=',num2str(alpherr),', abs error=',num2str(abs(zz(1,end)))])
        if (alpherr<tol)
            % right alph found, end search
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
        disp(['alph for current omg is: ', num2str(alph)])
        disp(['# of trials to find alph for this omg= ', num2str(itr)]); disp(' ')
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



