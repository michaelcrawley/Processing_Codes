function [data] = xLST(m,nmode,tr,omega,ialpha,x,y,U, varargin)
    %This code performs spatial linear stability theory calculations on PIV
    %data throughout the entire domain.  A curvefitting routine is used on
    %the PIV data to calculate R/theta at each streamwise location, which
    %is then passed to the LST function.
    %Version xLST 1.0
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
    %           use '-plot' to plot the scaled local eigenvalue along the
    %           streamwise axis
    %           use '-tag' to add tag to end of file name (ex. '-tag test')
    %           use '-tanh' for hyperbolic tangent velocity profile 
    %           use '-gaussian' for half-gaussian velocity profile

    multiWaitbar('Location:',0,'Color',[0.1 0.5 0.8]);%Initialize waitbars
    multiWaitbar('Iteration:',0,'Color',[0.2 0.9 0.3]); 
    multiWaitbar('Integration:',0,'Color',[0.8 0.4 0.9]); 
    
%  Data Formatting 
    if ~isempty(varargin) 
        if any(strcmp(varargin,'-standard') | strcmp(varargin,'-s'))
            x = x(:,1)';
            y = y(1,1:floor(end/2));
            U = U(:,1:floor(end/2))';
            Ureserve = U;
            Umax = max(U);
            for i = 1:length(x)
               U(:,i) = U(:,i)/Umax(i); 
            end       
        end
        if any(strcmp(varargin,'-tanh'))
            Ufun = @(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x))); %Velocity Profile
            newU1fun = @(x,o) 0.5*(1-(tanh(0.25*o*(1./x-x))).^2)*0.25*o.*(-1./x.^2-1); %First Derivative of Velocity Profile (with respect to x)
        elseif any(strcmp(varargin,'-gaussian'))
            Ufun = @(x,tau,A)exp(-tau.*(x-A).^2).*(x>A)+1.*(x<=A);
            newU1fun = @(x,o)(o(1).*(2*o(2) - 2*x))./exp(o(1).*(o(2) - x).^2);
        end
    else
        [M N] = size(U);
        if (length(x) ~= M && length(x) ~= N) || (length(y) ~= M && length(y) ~= N)
            error('Input matrix size mismatch');
        elseif (length(x) ~= N && length(x) == M)
            U = U';
        end        
    end
    
    if ~exist('Ufun','var')|| ~exist('newU1fun','var')
        Ufun = {@(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x)));@(x,b) (1+b*(x.^2)).^-2};
        newU1fun = {@(x,o) 0.5*(1-(tanh(0.25*o*(1./x-x))).^2)*0.25*o.*(-1./x.^2-1); @(x,o) -(4*o.*x)./(o*x.^2 + 1).^3};
    end
    
    data = EigenvalueCalc(m,nmode,tr,omega,ialpha,x,y,U,Ufun,newU1fun); 
    multiWaitbar('CLOSEALL');
    
    %  Optional Commmands
    if ~isempty(varargin)
        if any(strcmp(varargin,'-plot'))
            h(1) = figure;plot(x(1:length(data.Eigenvalue)),-imag(data.Eigenvalue)./data.Rt);title(['Streamwise Growth Rate: Mach(',num2str(m,2),') Tr(',num2str(tr,2),') AzMode(',num2str(nmode),') Omega(',num2str(omega),')']);xlabel('x/D');ylabel('-\alpha_i\theta');
            h(2) = figure;plot(x(1:length(data.Eigenvalue)),omega./real(data.Eigenvalue));title(['Disturbance Phase Speed: Mach(',num2str(m,2),') Tr(',num2str(tr,2),') AzMode(',num2str(nmode),') Omega(',num2str(omega),')']);xlabel('x/D');ylabel('\omega/\alpha_r');
        end
        if any(strcmp(varargin,'-save'))
            i = 0;
            dirname = ['Mach(',num2str(m,2),') Tr(',num2str(tr,2),') AzMode(',num2str(nmode),') Omega(',num2str(omega),')'];
            I = strmatch('-tag',varargin);
            if ~isempty(I)
                dirname = strcat(dirname,varargin{I}(5:end));
            end
            filename = strcat(dirname,'.mat');
            if exist(dirname,'file') ~= 7
                mkdir(pwd,dirname);
            end
            cd(dirname);
            filecheck = exist(filename,'file');
            while filecheck ~= 0        
                i = i+1;                
                if filecheck > 0
                    if i == 1
                        filename = filename(1:end-4);
                    else
                        filename = filename(1:end-7);
                    end
                    filename = strcat(filename,'(',num2str(i),').mat');
                end 
                filecheck = exist(filename,'file');
            end            
            if exist('h','var')
                for ii = 1:length(h)
                    saveas(h(ii),[filename(1:end-4),'_f',num2str(ii),'.fig']);
                end
                close(h);
            end
            if any(strcmp(varargin,'-perturb'))                
                up = perturb(data.Omega,data.Azimuthal_Mode,data.x,data.y,Ureserve,Umax,data.Rhalf,y,data.Eigenvalue,data.Eigenfunction,[filename(1:end-4) 'ptb.gif']);
            end
            tempstr = ['M',num2str(round(m*100)),'_Tr',num2str(round(tr*100)),'_AZ',num2str(nmode),'_Omg',num2str(round(omega*100)),'_',num2str(i)];
            tempvar = genvarname(tempstr);            
            eval([tempvar '=data;']);
            save(filename,tempstr);
            cd ..;
        end
    end        
end

function [data] = EigenvalueCalc(m,nmode,tr,omega,ialpha,x,y,U,Ufun, newU1fun)
%  Variable Initialization
    Tfun = @(x,U,tr) tr+(1-tr+((1.4-1)*m*m/2)).*(U)-((1.4-1)*m*m/2).*(U).^2; %Temperature Profile
    version = ['xLST 1.0';'LST 1.0 ']; 
   
       
%  Streamwise computation
    lalpha = zeros(1,length(x));
    alpherror = zeros(1,length(x));
    lphi = cell(1,length(x));
    
    lalpha(1) = ialpha;
    Rhalf = calcRhalf(x,y,U);
    
    for i = 1:length(x)
        multiWaitbar('Iteration:','Reset');
        multiWaitbar('Integration:','Reset');
        xt = y/Rhalf(i);
        yt = U(:,i);
        [Rt(i),~,~,~] = GeneralFit(xt,yt,@(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x))));
        Rt(i) = round(Rt(i)*1000)/1000;
        if ~iscell(Ufun)            
            [params{i},~,newUfun,rmserror(i)] = GeneralFit(xt,yt,Ufun);
            params{i} = round(params{i}*1E2)/1E2;
            cnewUfun{i} = newUfun;
            cnewU1fun{i} = newU1fun;
        else
            for ii = 1:length(Ufun)
                [estparams{ii}, ~,estnewUfun{ii},sse(ii)] = GeneralFit(xt,yt,Ufun{ii});
            end
            [~,I] = min(sse);
            params{i} = round(estparams{I}*1E2)/1E2;
            cnewUfun{i} = estnewUfun{I};
            rmserror(i) = sse(I);
            cnewU1fun{i} = newU1fun{I};
        end
        
        if i > 2
            lalpha(i) = interp1(x(1:i-1),lalpha(1:i-1),x(i),'spline'); 
        end        
        
        [lalpha(i) alpherror(i) lphi{i}] = LST(m,nmode,omega,params{i},tr,lalpha(find(lalpha,1,'last')),cnewUfun{i},cnewU1fun{i},Tfun);
        
        if isnan(lalpha(i))
            multiWaitbar('Location:','Color',[0.8 0.0 0.1]);
            multiWaitbar('Iteration:','Color',[0.8 0.0 0.1]);
            multiWaitbar('Iteration:','Value',1);
            multiWaitbar('Integration:','Color',[0.8 0.0 0.1]);
            multiWaitbar('Integration:','Value',1);
            pause(2);
            lalpha = lalpha(1,1:i);
            alpherror = alpherror(1,1:i);
            lphi = lphi(1,1:i);           
           break;
%         elseif imag(lalpha(i)) > 0
%             multiWaitbar('Location:','Color',[0.8 0.0 0.1]);
%             multiWaitbar('Iteration:','Color',[0.8 0.0 0.1]);
%             multiWaitbar('Iteration:','Value',1);
%             multiWaitbar('Integration:','Color',[0.8 0.0 0.1]);
%             multiWaitbar('Integration:','Value',1);
%             pause(2);
%            break;
         else
            multiWaitbar('Location:','Increment',1/length(x));
        end
    end
    date_run = date; 
    data = struct('Mach_Number',m,'Azimuthal_Mode',nmode,'Temp_Ratio',tr,'Omega',omega,'Rt',Rt,'Rhalf',Rhalf,'estparams',{params},'Eigenvalue',lalpha,'Eigenfunction',{lphi},'Code_Version',version,'Date_run',date_run,'Velocity_Profile',{cnewUfun},'x',x,'y',y,'U',U,'alpherror',alpherror,'RMSerror',rmserror); 
end

function [alph alpherror phi] = LST(m,nmode,omega,Rt,tr,alph,Ufun,U1fun,Tfun)
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
    tol = 5E-4; %tolerance for alpha
    y1 = .005;%Define range for velocity profile
    y2 = 4;
    u = Ufun(y1,Rt);
    t = Tfun(y1,u,tr);
    alpherror = 1;
    counter = 1;
    itrmax = 20;
    
    while (alpherror > tol) && (counter <= itrmax)
        clear phi;
        multiWaitbar('Iteration:','Increment',1/itrmax);
        %Initial values for integration (along centerline of jet)
        c = omega/alph;      
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
        phi(1,:) = z(1,:);
        phi(2,:) = y;

        %Calculate Integration Errors
        v = Ufun(y2,Rt);
        vmc = v-c;
        t2 = Tfun(y2,v,tr);
        wm=sqrt(alph*alph*(1-m*m*vmc*vmc/t2));

        er=wm*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-z(2,end)*besselh(nmode,wm*y2);
        dwma=((t2-m*m*vmc*vmc)*alph-m*m*omega*vmc)/(wm*t2);
        der=dwma*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+wm*z(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+...
            wm*z(1,end)*y2*dwma*(besselh(nmode-2,wm*y2)-2*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4-z(4,end)*besselh(nmode,wm*y2)-z(2,end)*dwma*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2;
        dalpha=-er/der;
        alpherror=abs(dalpha)/abs(alph);
        
        if (alpherror > tol) || isnan(alpherror)
            alph=alph+dalpha;
            if counter == itrmax || isnan(alpherror)
               alph = NaN;
               counter = itrmax+1;
            end
        end
        counter = counter+1;
    end
    multiWaitbar('Iteration:','Value',1);
    pause(1/60);
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

function [z,y]=odeRKF45(fun,z1,ystart,yend,abserr,inputs)
%   this code integrates a system of first order ordinary diferential
%   equations by runge-kutta-fehberg-45 method with automatic estimation
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

    if yend < ystart
        dx = -1;
    else
        dx = 1;
    end
    % step size
    h= 0.01*dx;      % h can be either positive or negative
    hmin=0.0001*dx;   % hmin must be positive
    hmax=0.02*dx;    % hmax must be positive

    % control parameters to avoid dead loop
    max1=200000;    % maximun number of valid integrations
    max2=400000;    % maximun number of all integrations
    i=0;          % number of all integration (including not valid integrations)
    j=1;          % number of valid integrations

    % results to be output
    y=[ ];
    y=[y,ystart];
    z=[ ];
    z=[z, z1];

    while abs(y(j)-yend) > eps
        multiWaitbar('Integration:','Value',abs((y(j)-ystart)/(yend-ystart)));
        if abs(y(j)-yend) < abs(h)
            h=abs(yend-y(j))*dx;
        end
        k1=h*feval(fun,y(j),z(:,j), inputs);
        k2=h*feval(fun,y(j)+a2*h, z(:,j)+b2*k1, inputs);
        k3=h*feval(fun,y(j)+a3*h, z(:,j)+b3*k1+c3*k2, inputs);
        k4=h*feval(fun,y(j)+a4*h, z(:,j)+b4*k1+c4*k2+d4*k3, inputs);
        k5=h*feval(fun,y(j)+a5*h, z(:,j)+b5*k1+c5*k2+d5*k3+e5*k4, inputs);
        k6=h*feval(fun,y(j)+a6*h, z(:,j)+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5, inputs);
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
            h=dx*hmax*h/abs(h);
        elseif abs(h)<hmin
            h=dx*hmin*h/abs(h);
        end
        % reach maximum number of integration, break
        i=i+1;
        if j==max1 || i==max2
            %disp(' ');
            %disp('error in RKF45: too many integrations needed');
            return
        end

    end
end


function [up u] = perturb(omega,nmode,x,y,um,Umax,Rhalf,yt,alpha,Eigenfunction,filename)
    x = x(1:length(alpha));
    um = um(:,1:length(alpha));
    Umax = Umax(1:length(alpha));
    A = 0.04*max(Umax); %set initial amplitude for disturbance

    t = (0:pi/32:2*pi)/omega;
    theta = [0 pi];
    
    %zmap = permute(repmat(x',[1 length(y) length(theta) length(t)]), [ 3 2 1 4]);
    thetamap = repmat(theta',[1 length(y) length(x) length(t)]);
    %rmap = repmat(y(end:-1:1), [length(theta) 1 length(x) length(t)]);
    tmap = permute(repmat(t,[length(theta) 1 length(x) length(y)]),[1 4 3 2]);
    alphamap = permute(repmat(cumtrapz(x,alpha).',[1 length(y) length(theta) length(t)]), [3 2 1 4]);   
    for i = 1:length(alpha)
        phi(:,i) = interp1(Rhalf(i)*Eigenfunction{i}(2,:),Eigenfunction{i}(1,:),yt,'spline');
    end   
    phimap = permute(repmat(phi,[1  1 length(theta) length(t)]),[3 1 2 4]);
    
    up = A*real(phimap.*exp(1i*(alphamap+nmode*thetamap-omega*tmap)));    
    up = [squeeze(up(1,:,:,:)); squeeze(up(2,end:-1:1,:,:))];
    um = repmat([um; um(end:-1:1,:)],[1 1 length(t)]);
       
    u = um+up;    
    
    save([filename(1:end-3) 'mat'],'up','u','phimap','alphamap');
    
    for i = 1:length(t)
        h(1) = figure;
        pcolor(x,[y -y(end:-1:1)],up(:,:,i)/max(Umax));shading interp;title(filename(1:end-4));xlabel('x/D');ylabel('y/D');colorbar;caxis([min(min(min(up)))/max(Umax) max(max(max(up)))/max(Umax)]);
        gif_add_frame(h(1),filename,10);
        close(h);
    end
end