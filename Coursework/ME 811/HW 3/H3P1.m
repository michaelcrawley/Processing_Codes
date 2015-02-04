function [phia] = H3P1(N,Pe)
    %Code completed by Michael Crawley for ME 811 HW#3

    x = 1/(2*N):1/N:1-1/(2*N);
    dx = mean(diff(x));
    
    bound0 = 0;
    bound1 = 1;
    
    gamma = 1;
    rhou = Pe*gamma;
    
    %%Analytic
    C = (exp(Pe)-1)^-1;
    phia = C*(exp(Pe*x)-1);
    
    %%First order upwind
    X0 = repmat(rhou+2*gamma/dx,N,1);
    X0([1 end]) = rhou+4*gamma/dx; %apply boundary conditions
    Xn1 = repmat(-0.5*(abs(rhou)+rhou)-gamma/dx,N,1);
    X1 = repmat(-0.5*(abs(rhou)-rhou)-gamma/dx,N,1);
    % i = 1 boundary conditions
    Xn1(1) = 0;
    X1(1) = X1(1)-gamma/dx/3;
    % i = N boundary conditions
    X1(end) = 0; 
    Xn1(end) = Xn1(end)-gamma/dx/3;
    X0(end) = 0.5*(abs(rhou)-rhou)+4*gamma/dx;
    
    
    X1 = circshift(X1,1);
    Xn1 = circshift(Xn1,-1);
    
    A = spdiags([Xn1 X0 X1], -1:1,N,N);
    b = zeros(N,1);
    b(1) = (0.5*(abs(rhou)+rhou)+8*gamma/dx/3)*bound0; %apply boundary condition at x = 0
    b(end) = (0.5*(abs(rhou)-rhou)-0.5*(abs(rhou)+rhou)+8*gamma/dx/3)*bound1; %apply boundary condition at x = 1
    phiFUD = TDMsolver(A,b);
    
    %%Second order upwind 
    b = zeros(N,1);
    X0 = repmat(0.75*(abs(rhou)+rhou)+0.75*(abs(rhou)-rhou)+2*gamma/dx,N,1);
    Xn1 = repmat(-(0.25*(abs(rhou)+rhou)+0.75*(abs(rhou)+rhou)+gamma/dx),N,1);
    X1 = repmat(-(0.25*(abs(rhou)-rhou)+0.75*(abs(rhou)-rhou)+gamma/dx),N,1);
    Xn2 = repmat(0.25*(abs(rhou)+rhou),N,1);
    X2 = repmat(0.25*(abs(rhou)-rhou),N,1);
    %i = 1 boundary conditions
    X0(1) = abs(rhou)+rhou+4*gamma/dx;
    Xn1(1) = 0;
    X1(1) = -(0.75*(abs(rhou)-rhou)+4*gamma/dx/3);
    Xn2(1) = 0;
    b(1) = -(abs(rhou)+rhou+0.5*(abs(rhou)-rhou)+8*gamma/dx/3)*bound0;
    %i = 2 boundary conditions
    Xn1(2) = -(1.25*(abs(rhou)+rhou)+gamma/dx);
    Xn2(2) = 0;
    b(2) = -0.5*(abs(rhou)+rhou)*bound0;
    %i = N boundary conditions
    X0(end) = (abs(rhou)-rhou)+4*gamma/dx;
    Xn1(end) = -(0.75*(abs(rhou)+rhou)+4*gamma/dx/3);
    X1(end) = 0;
    X2(end) = 0;
    b(end) = (abs(rhou)-rhou-0.5*(abs(rhou)+rhou)+8*gamma/dx/3)*bound1;
    %i = N-1 boundary conditions
    X1(end-1) = -(1.25*(abs(rhou)-rhou)+gamma/dx);
    X2(end-1) = 0;
    b(end-1) = -0.5*(abs(rhou)-rhou)*bound1;
    
    X1 = circshift(X1,1);
    X2 = circshift(X2,2);
    Xn1 = circshift(Xn1,-1);
    Xn2 = circshift(Xn2,-2);
    
    A = spdiags([Xn2 Xn1 X0 X1 X2],-2:2,N,N);
    phiSUD = A\b;    
    
    %%QUICK
    b = zeros(N,1);
    X0 = repmat((3/16)*(abs(rhou)+rhou)+(3/16)*(abs(rhou)-rhou)+2*gamma/dx,N,1);
    Xn1 = repmat(-((1/16)*(abs(rhou)+rhou)+(3/8)*(abs(rhou)+rhou)-(3/16)*(abs(rhou)-rhou)+gamma/dx),N,1);
    X1 = repmat((3/16)*(abs(rhou)+rhou)-(3/8)*(abs(rhou)-rhou)-(1/16)*(abs(rhou)-rhou)-gamma/dx,N,1);
    Xn2 = repmat((1/16)*(abs(rhou)+rhou),N,1);
    X2 = repmat((1/16)*(abs(rhou)-rhou),N,1);
    %i = 1 boundary conditions
    X0(1) = 0.5*(abs(rhou)+rhou)-(3/16)*(abs(rhou)-rhou)+4*gamma/dx;
    Xn1(1) = 0;
    X1(1) = (1/6)*(abs(rhou)+rhou)-(3/8)*(abs(rhou)-rhou)-4*gamma/dx/3;
    Xn2(1) = 0;
    b(1) = ((1/6)*(abs(rhou)+rhou)+rhou+8*gamma/dx/3)*bound0;
    %i = 2 boundary conditions
    X0(2) = (3/8)*(abs(rhou)+rhou)-(3/16)*(abs(rhou)-rhou)-(1/6)*(abs(rhou)+rhou)+(3/8)*(abs(rhou)-rhou)+2*gamma/dx;
    Xn1(2) = -(1/16)*(abs(rhou)+rhou)-0.5*(abs(rhou)+rhou)+(3/16)*(abs(rhou)-rhou)-gamma/dx;
    Xn2(2) = 0;
    b(2) = -((1/6)*(abs(rhou)+rhou))*bound0;
    %i = N boundary conditions
    X0(end) = 0.5*(abs(rhou)-rhou)-(3/16)*(abs(rhou)+rhou)+4*gamma/dx;
    Xn1(end) = (1/6)*(abs(rhou)-rhou)-(3/8)*(abs(rhou)+rhou)-4*gamma/dx/3;
    X1(end) = 0;
    X2(end) = 0;
    b(end) = ((1/6)*(abs(rhou)-rhou)-rhou+8*gamma/dx/3)*bound1;
    %i = N-1 boundary conditions
    X0(end-1) = (3/8)*(abs(rhou)+rhou)-(1/6)*(abs(rhou)-rhou)-(3/16)*(abs(rhou)+rhou)+(3/8)*(abs(rhou)-rhou)+2*gamma/dx;
    X1(end-1) = -(1/16)*(abs(rhou)-rhou)-0.5*(abs(rhou)-rhou)+(3/16)*(abs(rhou)+rhou)-gamma/dx;
    X2(end-1) = 0;
    b(end-1) = -((1/6)*(abs(rhou)-rhou))*bound1;
    
    X1 = circshift(X1,1);
    X2 = circshift(X2,2);
    Xn1 = circshift(Xn1,-1);
    Xn2 = circshift(Xn2,-2);
    
    A = spdiags([Xn2 Xn1 X0 X1 X2],-2:2,N,N);
    phiQUICK = A\b;
    
    
    %%Exponential
    Pel = Pe*dx;
    fp = Pel/(exp(Pel)-1);
    fm = Pel*exp(Pel)/(exp(Pel)-1);
        
    b = zeros(N,1);
    X0 = repmat(gamma*fm/dx+gamma*fp/dx,N,1);
    Xn1 = repmat(-gamma*fm/dx,N,1);
    X1 = repmat(-gamma*fp/dx,N,1);
    %i = 1 boundary conditions
    X0(1) = gamma*fm/dx+gamma*Pel/dx/(exp(Pel/2)-1);
    Xn1(1) = 0;
    b(1) = (rhou+gamma*Pel/dx/(exp(Pel/2)-1))*bound0;
    %i = N boundary conditions
    X0(end) = gamma*Pel*exp(Pel/2)/dx/(exp(Pel/2)-1)+gamma*fp/dx;
    X1(end) = 0;
    b(end) = (-rhou+gamma*Pel*exp(Pel/2)/dx/(exp(Pel/2)-1))*bound1;
    
    X1 = circshift(X1,1);
    Xn1 = circshift(Xn1,-1);
    
    A = spdiags([Xn1 X0 X1], -1:1,N,N);
    phiEXP = TDMsolver(A,b);
    
    %%Plot Results
    h(1) = figure;
    plot(x,phia);xlabel('x');ylabel('\Phi_a');title(['Analytic Solution for Pe = ',num2str(Pe)]);
    h(2) = figure;
    plot(x,phiFUD'-phia,'k',x,phiSUD'-phia,'-sk',x,phiQUICK'-phia,'--k',x,phiEXP'-phia,'-.k');
    legend('1^s^t Order Upwind','2^n^d Order Upwind', 'QUICK','Exponential','Location','Best');
    xlabel('x');
    ylabel('\Phi_n-\Phi_a');
    title(['Error in numerical solution for N = ',num2str(N),', Pe = ', num2str(Pe)]);
    saveas(h(2),['N',num2str(N),' Pe',strrep(num2str(Pe),'.','_')],'fig');
    saveas(h(2),['N',num2str(N),' Pe',strrep(num2str(Pe),'.','_')],'png');
    saveas(h(1),['Phia Pe',strrep(num2str(Pe),'.','_')],'fig');
    saveas(h(1),['Phia Pe',strrep(num2str(Pe),'.','_')],'png');
    close(h);
end