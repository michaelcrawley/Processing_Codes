function [lalpha lphi] = xLSTv3(m,nmode,tr,omega,ialpha,x,y,U, varargin)
    %This code performs spatial linear stability theory calculations on PIV
    %data throughout the entire domain.  A curvefitting routine is used on
    %the PIV data to calculate R/theta at each streamwise location, which
    %is then passed to the LST function.
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
    %           use '-tanh' for hyperbolic tangent profile (default)
    %           use '-gaussian' for one sided gaussian profile
    %           use '-save' to save all variables to a .mat file.  
    
    %need to add save function, plot function
              
    %Data formatting & Curve-fitting
    if ~isempty(varargin) 
        if any(strcmp(varargin,'-standard') | strcmp(varargin,'-s'))
            x = x(3:127,1);
            y = y(1,1:44);
            U = U(3:127,1:44)';
            Umax = max(U);
            for ii = 1:length(x)
               U(:,ii) = U(:,ii)/Umax(ii); 
            end 
        end
        if any(strcmp(varargin,'-tanh')) 
            Ufun = @(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x)));
            U1fun = @(x,Rt) 0.125*Rt.*(tanh(0.25*Rt.*(x - 1./x)).^2 - 1).*(1./x.^2 + 1);
        elseif any(strcmp(varargin,'-gaussian'))
            Ufun = @(x,tau,A)exp(-tau.*(x-A).^2).*(x>A)+1.*(x<=A);
            U1fun = @(x,tau,A)(tau.*(2*A - 2*x))./exp(tau.*(A - x).^2);
        end
        
        for i = 1:length(varargin)
            if isa(varargin{i},'function_handle')
                Ufun = varargin{i};
            end
        end
    elseif isempty(varargin) || ~any(strcmp(varargin,'-standard') | strcmp(varargin,'-s'))
        [M N] = size(U);
        if (length(x) ~= M && length(x) ~= N) || (length(y) ~= M && length(y) ~= N)
            error('Input matrix size mismatch');
        elseif (length(x) ~= N && length(x) == M)
            U = U';
        end
    end    
    if ~exist('Ufun','var')%use default (tanh) velocity profile       
        warning('Curve-fit for velocity profile not specified, using default tanh profile'); %#ok<WNTAG>
        Ufun = @(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x)));
        U1fun = @(x,Rt) 0.125*Rt.*(tanh(0.25*Rt.*(x - 1./x)).^2 - 1).*(1./x.^2 + 1);
    end
    
    %Initialize variables
    lalpha = zeros(1,length(x));
    lphi = cell(1,length(x));
    lalpha(1) = ialpha;
    
    %Streamwise computation    
    Rhalf = calcRhalf(x,y,U);
    for i = 1:length(x)
        xt = y/Rhalf(i);
        yt = U(:,i);
        y1t = NumericalDerivative(1,4,mean(diff(xt)),yt);
        [Uparams,~, newUfun] = GeneralFit(xt,yt,Ufun);
        [~, ~, newU1fun] = GeneralFit(xt,y1t,U1fun);
        if i > 2
            lalpha(i) = interp1(x(1:i-1),lalpha(1:i-1),x(i),'spline'); 
        end
        [lalpha(i),lphi{i}] = LSTv2(m,nmode,tr,newUfun,Uparams,newU1fun,Uparams,xt,omega,lalpha(end));
        if isnan(lalpha(i))
           warning(['Local Eigenvalue not found at x/D = ',num2str(x(i)),'; ending program']); %#ok<WNTAG>
           break;
        elseif imag(lalpha(i)) > 0
           warning(['Local Eigenvalue at x/D = ',num2str(x(i)),', is damped; ending program']); %#ok<WNTAG>
           break;
        else
            disp(['Processing completed for x/D = ',num2str(x(i)),', found eigenvalue = ',num2str(lalpha(i))]);
        end
    end    
    phasevelocity = real(omega./lalpha);
    
    if ~isempty(varargin)
       if any(strcmp(varargin,'-plot'))
           %add plot commands
           figure;plot(x(1:i),phasevelocity);xlabel('x/D');ylabel('c_r');title('Phase Velocity');
       end
       if any(strcmp(varargin,'-save'))
            i = 1;
            filecheck = 1;
            filename = strcat(pwd,'\Mach(',num2str(m,2),') AzMode(',num2str(nmode),') Tr(',num2str(tr,2),').mat');
            while filecheck ~= 0        
                filecheck = exist(filename,'file');
                if filecheck > 0
                    filename = strcat(pwd,'\Mach(',num2str(m,2),') AzMode(',num2str(nmode),') Tr(',num2str(tr,2),')(',num2str(i),').mat');
                end
                i = i+1;
            end
            clear i ii xt yt Ufit M N filecheck varargin
            date_run = date; %#ok<NASGU>
            mfilename = 'xLSTv3'; %#ok<NASGU>
            save(filename);
       end       
    end
end