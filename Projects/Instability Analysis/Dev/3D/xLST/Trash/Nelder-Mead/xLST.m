function[] = xLST(m,nmode,tr,omega,ialpha,x,y,U, varargin)
    %This code performs spatial linear stability theory calculations on PIV
    %data throughout the entire domain.  A curvefitting routine is used on
    %the PIV data to calculate R/theta at each streamwise location, which
    %is then passed to the LST function.
    %Version xLST NM.1.0.0.0
    %Inputs:
    %m : Mach Number
    %nmode : Azimuthal Mode
    %tr : Temperature Ratio
    %omega : perturbation angular frequency
    %ialpha : initial guess for the LST eigenvalue
    %x : streamwise locations
    %y : cross streamwise locations
    %U : mean velocity profile
    %options :  use '-standard' or '-s' for standard PIV data input
    %           use '-save' to save data to a .mat file
    
%  Data Formatting 
    if ~isempty(varargin) 
        if any(strcmp(varargin,'-standard') | strcmp(varargin,'-s'))
            x = x(3:127,1);
            y = y(1,1:44);
            U = U(3:127,1:44)';
            Umax = max(U);
            for i = 1:length(x)
               U(:,i) = U(:,i)/Umax(i); 
            end       
        end
    else
        [M N] = size(U);
        if (length(x) ~= M && length(x) ~= N) || (length(y) ~= M && length(y) ~= N)
            error('Input matrix size mismatch');
        elseif (length(x) ~= N && length(x) == M)
            U = U';
        end
    end
    data = EigenvalueCalc(m,nmode,tr,omega,ialpha,x,y,U);  
    %  Optional Commmands
    if ~isempty(varargin)
        if any(strcmp(varargin,'-plot'))
            figure;plot(x,-imag(data.Eigenvalue));title('Streamwise Growth Rate');xlabel('x/D');ylabel('-\alpha_i');
        end
        if any(strcmp(varargin,'-save'))
            i = 1;
            filecheck = 1;
            filename = strcat(pwd,'\Mach(',num2str(m,2),') Tr(',num2str(tr,2),') AzMode(',num2str(nmode),') Omega(',num2str(omega),').mat');
            while filecheck ~= 0        
                filecheck = exist(filename,'file');
                if filecheck > 0
                    filename = strcat(pwd,'\Mach(',num2str(m,2),') Tr(',num2str(tr,2),') AzMode(',num2str(nmode),') Omega(',num2str(omega),')(',num2str(i),').mat');
                end
                i = i+1;
            end
            tempstr = ['M',num2str(round(m*100)),'_Tr',num2str(round(tr*100)),'_AZ',num2str(nmode),'_Omg',num2str(round(omega*100))];
            tempvar = genvarname(tempstr);            
            eval([tempvar '=data;']);
            save(filename,tempstr);
        end     
    end
end

function [data] = EigenvalueCalc(m,nmode,tr,omega,ialpha,x,y,U)
%  Variable Initialization
    Ufun = @(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x))); %Velocity Profile
    U1fun = @(x,Rt) 0.5*(1-(tanh(0.25*Rt*(1./x-x))).^2)*0.25*Rt.*(-1./x.^2-1); %First Derivative of Velocity Profile (with respect to x)
    Tfun = @(x,U,tr) tr+(1-tr+((1.4-1)*m*m/2)).*(U)-((1.4-1)*m*m/2).*(U).^2; %Temperature Profile
    Rt = zeros(1,length(x));
    version = 'xLST NM.1.0.0.0'; 
    h = waitbar(0,'Starting Processing...');
       
%  Streamwise computation
    lalpha(1) = ialpha;
    Rhalf = calcRhalf(x,y,U);
    for i = 1:length(x)
        xt = y/Rhalf(i);
        yt = U(:,i);
        [Rt(i),~,~,rmserror(i)] = GeneralFit(xt,yt,Ufun);
        if i > 2
            lalpha(i) = interp1(x(1:i-1),lalpha(1:i-1),x(i),'spline'); 
        end
        lalpha(i)= fzero(@errorcalc,lalpha(end));
        [alpherror(i),lphi(i,:)] = errorcalc(lalpha(i));
        if isnan(lalpha(i))
            waitbar(1,h,['Local Eigenvalue not found at x/D = ',num2str(x(i)),'; ending program']);
            pause(5);
           break;
        elseif imag(lalpha(i)) > 0
            waitbar(1,h,['Local Eigenvalue at x/D = ',num2str(x(i)),', is damped; ending program']);
            pause(5);
           break;
        else
            waitbar(i/length(x),h,['Eigenvalue found @ x/D = ', num2str(round(10*x(i))/10)]);
        end
    end
    function[error, phi] = errorcalc(alph)
        %Calculate Integration Errors
        [z inty] = LST(m,nmode,omega,Rt(i),tr,alph,Ufun,U1fun,Tfun);
        y2 = inty(end);
        c = omega/alph;
        v = Ufun(y2,Rt(i));
        vmc = v-c;
        t2 = Tfun(y2,v,tr);
        wm=sqrt(alph*alph*(1-m*m*vmc*vmc/t2));
        
        er=wm*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-z(2,end)*besselh(nmode,wm*y2);
        dwma=((t2-m*m*vmc*vmc)*alph-m*m*omega*vmc)/(wm*t2);
        der=dwma*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+wm*z(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+...
            wm*z(1,end)*y2*dwma*(besselh(nmode-2,wm*y2)-2*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4-z(4,end)*besselh(nmode,wm*y2)-z(2,end)*dwma*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2;
        dalpha=-er/der;
        error=abs(dalpha)/abs(alph);
        phi = z(1,:);
    end
    date_run = date; 
    data = struct('Mach_Number',m,'Azimuthal_Mode',nmode,'Temp_Ratio',tr,'Omega',omega,'Rt',Rt,'Eigenvalue',lalpha,'Eigenfunction',{lphi},'Code_Version',version,'Date_run',date_run,'Velocity_Profile',Ufun,'x',x,'y',y,'U',U,'alpherror',alpherror,'RMSerror',rmserror); 
    close(h);
end

function [z y] = LST(m,nmode,omega,Rt,tr,alph,Ufun,U1fun,Tfun)
    %Inputs:
    %m: Mach number
    %nmode: Azimuthal mode
    %omgR: Omega range
    %alph: Initial guess for alpha
    %Rt: R/theta
    %tr: temperature ratio (T_infinity/T_jet)
    %Ufun: Velocity Profile
    %Tfun: Temperature Profile

    abserr=5.e-4;%tolerance for ode solver  
    y1 = .005;%Define range for velocity profile
    y2 = 4;

    c = omega/alph;

    %Initial values for integration (along centerline of jet)
    u = Ufun(y1,Rt);
    t = Tfun(y1,u,tr);           
    umc=u-c;
    wl=sqrt(alph*alph*(1.-m*m*umc*umc/t));
    dwla=((t-m*m*umc*umc)*alph-m*m*omega*umc)/(wl*t);            
    z0(1,1)= besselj(nmode,wl*y1);
    z0(2,1)= wl*(besselj(nmode-1,wl*y1)-besselj(nmode+1,wl*y1))/2;
    z0(3,1)= z0(2,1)*y1*dwla/wl;
    z0(4,1)= z0(3,1)/y1+wl*(besselj(nmode-2,wl*y1)-2.*besselj(nmode,wl*y1)+besselj(nmode+2,wl*y1))*y1*dwla/4;

    %Integration (centerline to inf)
    params = struct('alph',alph,'omega',omega,'m',m,'nmode',nmode,'Rt',Rt,'tr',tr,'Ufun',Ufun,'U1fun',U1fun);
    [z,y]=odeRKF45(@derive,z0,y1,y2,abserr,params); 
end

function dz=derive(y,z,params)
    %Inputs:
    % y : independent variable
    % z : dependent variable
    % params : additional parameters
    %Outputs:
    % dz = dz/dy at y
    
    %Parameter recoupling
    alpha = params.alph;
    omega = params.omega;
    m = params.m;
    nmode = params.nmode;
    Rt = params.Rt;
    tr = params.tr;
   
    c = omega/alpha;
    
    %Velocity Profile
    u=params.Ufun(y,Rt);
    u1=params.U1fun(y,Rt);
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