function [Rt lalpha lphi] = xLSTv1(m,nmode,tr,omega,ialpha,x,y,U, varargin)
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
    %           use '-save' to save all to a .mat file (not yet finished)
    %           use '-plot' to output predefined plots at conclusion of
    %           calculations (not yet finished)
    
    h = waitbar(0,'Starting Processing...');
    version = 'xLST 1.0.0.0';
    Ufun = @(x,Rt) 0.5*(1+tanh(0.25*Rt*(1./x-x)));
           
    %Data formatting
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
    
    %Streamwise computation
    lalpha(1) = ialpha;
    Rhalf = calcRhalf(x,y,U);
    for i = 1:length(x)
        xt = y/Rhalf(i);
        yt = U(:,i);
        Rt(i) = GeneralFit(xt,yt,Ufun);
        if i > 2
            lalpha(i) = interp1(x(1:i-1),lalpha(1:i-1),x(i),'spline'); 
        end
        [lalpha(i),lphi{i}] = LSTv1(m,nmode,omega,lalpha(end),Rt(i),tr);
        if isnan(lalpha(i))
%           warning(['Local Eigenvalue not found at x/D = ',num2str(x(i)),'; ending program']); %#ok<WNTAG>
            waitbar(1,h,['Local Eigenvalue not found at x/D = ',num2str(x(i)),'; ending program']);
            pause(5);
           break;
        elseif imag(lalpha(i)) > 0
%           warning(['Local Eigenvalue at x/D = ',num2str(x(i)),', is damped; ending program']); %#ok<WNTAG>
            waitbar(1,h,['Local Eigenvalue at x/D = ',num2str(x(i)),', is damped; ending program']);
            pause(5);
           break;
        else
%             disp(['Processing completed for x/D =',num2str(x(i)),', found eigenvalue = ',num2str(lalpha(i))]);
            waitbar(i/length(x),h,['Eigenvalue found @ x/D = ', num2str(x(i))]);
        end
    end
    
    if ~isempty(varargin)
        if any(strcmp(varargin,'-plot'))
            %add plot commands
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
            date_run = date; %#ok<NASGU>
            tempstr = ['M',num2str(round(m*10)),'_Tr',num2str(round(tr*100)),'_AZ',num2str(nmode),'_Omg',num2str(round(omega*10))];
            tempvar = genvarname(tempstr);
            data = struct('Mach_Number',m,'Azimuthal_Mode',nmode,'Temp_Ratio',tr,'Omega',omega,'Rt',Rt,'Eigenvalue',lalpha,'Eigenfunction',{lphi},'Code_Version',version,'Velocity_Profile',Ufun,'x',x,'y',y,'U',U); 
            eval([tempvar '=data;']);
            save(filename,tempstr);
        end     
    end
    close(h);
end