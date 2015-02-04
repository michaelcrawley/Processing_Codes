function [] = AxiJet(m,nmode,omgR,alph,Rt,tr)
    %Inputs:
    %m: Mach number
    %nmode: Azimuthal mode
    %omgR: Omega range
    %alph: Initial guess for alpha
    %Rt: R/theta
    %tr: temperature ratio (T_infinity/T_jet)
    direct = 'Z:\My Documents\Matlab\Instability Analysis\';
    version = 1.0;
    tic;
    %Velocity Profile
    U = @(y,Rt) 0.5*(1+tanh(0.25*Rt*(1./y-y)));
    %Temperature Profile
    T = @(y,Rt,tr) tr+(1-tr+((1.4-1)*m*m/2))*(U(y,Rt))-((1.4-1)*m*m/2)*(U(y,Rt))^2;
    
    y1 = .005;%Define range for velocity profile
    y2 = 6;    
    itrmax = 50; %Max iterations to find alpha
    tol=1.e-5; %tolerance for found alpha
    abserr=1.e-5;%tolerance for ode solver  
    
    l = length(omgR);
    
    multiWaitbar('Frequency Range:',0,'Color',[0.1 0.5 0.8]);%Initialize waitbars
    multiWaitbar('Iteration:',0,'Color',[0.2 0.9 0.3]);  
    
    for i = 1:l
        omega = omgR(i);
        if i>2 %Use spline curve to guess next value for alpha
            alph = interp1(omgR(1:i-1),talpha(1:i-1),omega,'spline','extrap');
        end
        multiWaitbar('Iteration:','Reset');
        multiWaitbar('Iteration:','Color',[0.2 0.9 0.3]); %Reset Iteration waitbar for next omega
        for j = 1:itrmax
            multiWaitbar('Iteration:','Increment',1/itrmax);
            alpha{i}(j) = alph;
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
            v = U(y1,Rt);
            t2 = T(y1,Rt,tr);            
            vmc = v-c;
            wm=sqrt(alph*alph*(1-m*m*vmc*vmc/t2));
            
            %Integration (centerline to inf)
            [z,y{i}{j}]=odeRKF45(@derive,z0,y1,y2,abserr,[alph omega m nmode Rt tr]); 

            %Calculate Errors            
            er{i}(j)=wm*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-z(2,end)*besselh(nmode,wm*y2);
            dwma{i}(j)=((t2-m*m*vmc*vmc)*alph-m*m*omega*vmc)/(wm*t2);
            der{i}(j)=dwma{i}(j)*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+wm*z(3,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2+...
                wm*z(1,end)*y2*dwma{i}(j)*(besselh(nmode-2,wm*y2)-2*besselh(nmode,wm*y2)+besselh(nmode+2,wm*y2))/4-z(4,end)*besselh(nmode,wm*y2)-z(2,end)*dwma{i}(j)*y2*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2;
            dalpha{i}(j)=-er{i}(j)/der{i}(j);
            alphaerr{i}(j)=abs(dalpha{i}(j))/abs(alph);
            
            phi{i}{j} = z(1,:);
            psi{i}{j} = z(3,:);
           
            
            if j>itrmax/2
                multiWaitbar('Iteration:','Color',[0.8 0.0 0.1]); %Change color to red to signify difficulty converging on solution
            end          
            if (alphaerr{i}(j)<tol)% right alpha found, end search
                talpha(i) = alpha{i}(j);
                multiWaitbar('Iteration:','Value',1);
                pause(1/60);
                break            
            else% new alph for next trial
                flag = 1;
                while flag == 1
                    alphtemp = alph+dalpha{i}(j);
                    c = omega/alphtemp;
                    umc=u-c;
                    vmc = v-c;
                    wl=sqrt(alphtemp*alphtemp*(1.-m*m*umc*umc/t));
                    wm=sqrt(alphtemp*alphtemp*(1-m*m*vmc*vmc/t2));
                    dwla=((t-m*m*umc*umc)*alphtemp-m*m*omega*umc)/(wl*t);            
                    z0(1,1)= besselj(nmode,wl*y1);
                    z0(2,1)= wl*(besselj(nmode-1,wl*y1)-besselj(nmode+1,wl*y1))/2;
                    z0(3,1)= z0(2,1)*y1*dwla/wl;
                    z0(4,1)= z0(3,1)/y1+wl*(besselj(nmode-2,wl*y1)-2.*besselj(nmode,wl*y1)+besselj(nmode+2,wl*y1))*y1*dwla/4;
                    [z,ytemp]=odeRKF45(@derive,z0,y1,y2,abserr,[alphtemp omega m nmode Rt tr]);                 
                    
                    ertemp=wm*z(1,end)*(besselh(nmode-1,wm*y2)-besselh(nmode+1,wm*y2))/2.-z(2,end)*besselh(nmode,wm*y2);
                    
                    if abs(ertemp) > abs(er{i}(j))
                        dalpha{i}(j) = 0.5*dalpha{i}(j);
                    else
                        flag = 0;
                    end
                end
                alph=alph+dalpha{i}(j);
                if imag(alph) > 0.1
                    alph = complex(real(alph),-imag(alph));
                end
            end
        end
        if imag(alph) > 0
            multiWaitbar('Frequency Range:','Color',[0.8 0.0 0.1]);
        end
        if j == itrmax %End program if convergence on alpha is not met
            multiWaitbar('Frequency Range:','Color',[0.8 0.0 0.1]);
            multiWaitbar('Frequency Range:','Value',1);
            disp('search failed: No alph found for this omg');
            pause(1);
            break
        end
        multiWaitbar('Frequency Range:','Increment',1/l);
    end
    compute_time = toc/3600;
    data = struct('Mach_Number',m,'Azimuthal_Mode',nmode,'theta',1/Rt,'Temp_Ratio',tr,'Omega_Range',{omgR(1:length(talpha))},'True_Alpha',{talpha},'Alpha_Trace',{alpha},'y',{y},'Phi',{phi},'Psi',{psi},'Error',{er},'D_error',{der},'D_alpha',{dalpha},'Alpha_Error', {alphaerr}, 'Compute_time',compute_time,'Code_Version',version); 
    savedir = uigetdir(direct,'Specify Save Directory for Data');
    filename = strcat(savedir,'\Mach(',num2str(m,2),') AzMode(',num2str(nmode),') Rt(',num2str(Rt),') Tr(',num2str(tr,2),').mat');
    save(filename,'data');
    multiWaitbar('CLOSEALL');
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

