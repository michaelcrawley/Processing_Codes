function [master] = PlanarShearLayer(omgR,alph,L,W)
    l = length(omgR);

    y1 = -8;%Define range for velocity profile
    y2 = 8;    
    itrmax = 100; %Max iterations to find alpha
    tol=1.e-5; %tolerance for found alpha
    abserr=1.e-5;%tolerance for ode solver
    
    talpha = zeros(size(omgR));
    
    multiWaitbar('Frequency Range:',0,'Color',[0.1 0.5 0.8]);%Initialize waitbars
    multiWaitbar('Iteration:',0,'Color',[0.2 0.9 0.3]);
    
    for i = 1:l
        omega = omgR(i);
        if i>2 %Use spline curve to guess next value for alpha
            salphr=interp1(omgR(1:i-1),real(talpha(1:i-1)),omega,'spline','extrap');
            salphi=interp1(omgR(1:i-1),imag(talpha(1:i-1)),omega,'spline','extrap');
            alph=complex(salphr,salphi);
        end
        multiWaitbar('Iteration:','Reset');
        multiWaitbar('Iteration:','Color',[0.2 0.9 0.3]); %Reset Iteration waitbar for next omega
        for j = 1:itrmax
            multiWaitbar('Iteration:','Increment',1/itrmax);
            alpha{i}(j) = alph;
            c = omega/alph;
            z0=[exp(alph*y1); alph*exp(alph*y1); y1*exp(alph*y1); alph*y1*exp(alph*y1)+exp(alph*y1)];%Initial values for integration
            zend = [exp(-alph*y2); -alph*exp(-alph*y2); -y2*exp(-alph*y2); alph*y2*exp(-alph*y2)-exp(-alph*y2)];
            
            [z,y{i}{j}]=odeRKF45(@deriv,z0,y1,y2,abserr,[alph c L W]); %Integration (left to right)
            
            phi{i}{j} = z(1,:);
            psi{i}{j} = z(3,:);
            
            er{i}(j)=alph*z(1,end)+z(2,end);
            der{i}(j)=z(1,end)+alph*z(3,end)+z(4,end);
            dalpha{i}(j)=(-er{i}(j))/der{i}(j);
            alphaerr{i}(j)=abs(dalpha{i}(j))/abs(alpha{i}(j));
            
            if j>itrmax/2
                multiWaitbar('Iteration:','Color',[0.8 0.0 0.1]); %Change color to red to signify difficulty converging on solution
            end          
            if (alphaerr{i}(j)<tol)% right alpha found, end search
                talpha(i) = alpha{i}(j);
                multiWaitbar('Iteration:','Value',1);
                pause(1/60);
                break            
            else% new alph for next trial                
                alph=alph+dalpha{i}(j);
            end
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
    master = struct('Omega_Range',{omgR},'Alpha',{talpha},'Alpha_Trace',{alpha},'y',{y},'Phi',{phi},'Psi',{psi},'Error',{er},'D_error',{der},'D_alpha',{dalpha},'Alpha_Error', {alphaerr}); 
    multiWaitbar('CLOSEALL');
end

function dz=deriv(y,z,inputs)
    % y - independent variable
    % z - dependent variable
    % alpha - complex wave number
    % c = omega/alpha
    % dz = dz/dy at y
    alpha = inputs(1);
    c = inputs(2);
    L = inputs(3);
    W = inputs(4);
    
    u=1+L*tanh(y)-W*exp(-log(2)*y.^2);
    u2=-2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2);
    dz=zeros(4,1);
    cmc=u-c;
    dz(1)=z(2);                        
    dz(2)=(alpha*alpha+u2/cmc)*z(1) ;  
    dz(3)=z(4);						   
    dz(4)=(alpha*alpha+u2/cmc)*z(3)+z(1)*(2*alpha-c*u2/(cmc*cmc*alpha));
end

