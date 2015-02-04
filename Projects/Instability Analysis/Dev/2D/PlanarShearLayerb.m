function [master] = PlanarShearLayerb(omgR,alph,L,W)
    %This function uses a second order finite difference scheme, as opposed
    %to the Runge-Kutta method, to solve the ODE
    l = length(omgR);
    y1 = -10;%Define range for velocity profile
    y2 = 10;
    h = 0.01;
    y = y1:h:y2;
    ly = length(y);
    itrmax = 100; %Max iterations to find alpha
    tol=1.e-1; %tolerance for found alpha
    
    U = 1+L*tanh(y)-W*exp(-log(2)*y.^2);
    U2 = (-2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2));
    A = ones(1,ly);
    B = zeros(1,ly);
    D = zeros(1,ly);
    E = zeros(1,ly);
    F = ones(1,ly);
    G = zeros(1,ly);
    J = zeros(1,ly);
    f = zeros(1,ly);
    g = zeros(1,ly);
    
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
            
            sigma = [exp((y1-h)*alph) exp(-(y2+h)*alph)];%Boundary Conditions for Phi
            lambda = [(y1-h)*exp((y1-h)*alph) -(y2+h)*exp(-(y2+h)*alph)];%Boundary Conditions for Psi
            C = -(alph^2+U2./(U-omega/alph));
            H = C;
            K = -(2*alph-(omega*U2)./((U*alph-omega).^2));
            [phi{i}{j},psi{i}{j}] = SOFDM(A,B,C,D,E,F,G,H,J,K,f,g,sigma,lambda,h); %Integration based off of second order method
            phidy = NumericalDerivative(1,3,h,phi{i}{j}.');
            psidy = NumericalDerivative(1,3,h,psi{i}{j}.');
            
            er = mean([ -alph*phi{i}{j}(1)+phidy(1); alph*phi{i}{j}(end)+phidy(end)])
            der = mean([ -phi{i}{j}(1)-alph*psi{i}{j}(1)+psidy(1); phi{i}{j}(end)+alph*psi{i}{j}(end)+psidy(end)])
            dalpha = (-er)/der            
            alphaerr=abs(dalpha)/abs(alph)
            
            if j>itrmax/2
                multiWaitbar('Iteration:','Color',[0.8 0.0 0.1]); %Change color to red to signify difficulty converging on solution
            end          
            if (alphaerr<tol)% right alpha found, end search
                talpha(i) = alpha{i}(j);
                multiWaitbar('Iteration:','Value',1);
                pause(1/60);
                break            
            else% new alph for next trial                
                alph=alph+mean(dalpha)
            end
        end
        if j == itrmax %End program if convergence on alpha is not met
            multiWaitbar('Frequency Range:','Color',[0.8 0.0 0.1]);
            multiWaitbar('Frequency Range:','Value',1);
            disp('search failed: No alpha found for this omega');
            pause(1);
            break
        end
        multiWaitbar('Frequency Range:','Increment',1/l);
    end
    master = struct('Omega_Range',{omgR},'Alpha',{talpha},'Alpha_Trace',{alpha},'y',{y},'Phi',{phi},'Psi',{psi},'Error',{er},'D_error',{der},'D_alpha',{dalpha},'Alpha_Error', {alphaerr}); 
    multiWaitbar('CLOSEALL');
end