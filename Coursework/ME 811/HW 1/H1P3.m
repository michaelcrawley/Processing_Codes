function [Fluxl Fluxr Fluxnet x Error] = H1P3(N)
% Completed for HW1, Problem 3 for ME 811 by Michael Crawley
% Code performs Finite volume method to solve ODE
    

    %set boundary conditions
    boundl = 0;
    boundr = 1;    

    %generate mesh
    s = 1.02;    
    dx0 = (1-s)/(s-s^N);
    dx = (s.^(1:N-1))*dx0;
    x = [0 cumsum(dx)];
    x = (x(1:end-1)+x(2:end))/2;
    
    %calculate analytic solution
    phia = (1/3)*(x.^3)-(1/2)*(x.^2)+(7/6)*x;
    
    %generate matrices for linear equations
    coefs = [2./(dx(1:end-2)'+dx(2:end-1)') -(2./(dx(1:end-2)'+dx(2:end-1)')+2./(dx(2:end-1)'+dx(3:end)')) 2./(dx(2:end-1)'+dx(3:end)')];
    coefs = [0 -(2./(dx(1)+dx(2))+2/dx(1)) 2./(dx(1)+dx(2)); coefs; 2./(dx(end-1)+dx(end)) -(2./(dx(end-1)+dx(end))+2/dx(end)) 0];
    A = spdiags(fliplr(coefs),-1:1,N-1,N-1)';
    b = (2*x'-1).*dx';
    b(1) = b(1)-2*boundl/dx(1);
    b(end) = b(end)-2*boundr/dx(end);
    
    %calculate numeric solution
    phin = TDMsolver(A,b)';   
    
    %calculate fluxes at bounds and error in solution
    Fluxl = -2*(phin(1)-boundl)/dx(1);
    Fluxr = -2*(boundr-phin(end))/dx(end);
    Fluxgen = trapz([0 x 1],2*[0 x 1]-1);
    Fluxnet = Fluxr-Fluxl-Fluxgen;
    Error = abs(phia-phin);
    
    %plot
    plot(x,Error);title(['Error in FV solution, N = ',num2str(N)]);xlabel('x');ylabel('error (abs(analytic-numeric))');xlim([0 1]);

end