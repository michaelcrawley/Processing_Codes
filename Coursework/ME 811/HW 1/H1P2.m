function [Fluxl Fluxr Fluxnet x Error] = H1P2(N)
% Completed for HW1, Problem 2 for ME 811 by Michael Crawley
% Code performs Finite difference method to solve ODE
    
    %generate mesh
    s = 1.02;
    dx0 = (1-s)/(s-s^N);
    dx = (s.^(1:N-1))*dx0;
    x = [0 cumsum(dx)];
    
    %calculate analytical solution
    phia = (1/3)*(x.^3)-(1/2)*(x.^2)+(7/6)*x;
    
    %generate matrices for linear equations
    coefs = zeros(3,N-2);
    for j = 2:length(x)-1
        coefs(:,j-1) = TSE(1,1,x(j-1:j+1),2);
    end
    A = spdiags(flipud([[0 1 0]' coefs [0 1 0]'])',-1:1,N,N)';
    b = [0 2*x(2:end-1)-1 1]';
    
    %calculate numeric solution
    phin = TDMsolver(A,b)';

    %calculate fluxes and error in solution
    Fluxl = -phin(1:3)*TSE(0,2,x(1:3),1);
    Fluxr = -phin(end-2:end)*TSE(2,0,x(end-2:end),1);
    Fluxgen = trapz(x,2*x-1);
    Fluxnet = Fluxr-Fluxl-Fluxgen;
    Error =abs(phia-phin);
    
    %plot
    plot(x,Error);title(['Error in FD solution, N = ',num2str(N)]);xlabel('x');ylabel('error (abs(analytic-numeric))');xlim([0 1]);

end