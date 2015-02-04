function [Tn Ta] = P1(alpha, h, k, L, N, To, Tamb,tfinal)
%% Variable initialization
    x = 0:L/(N-1):L; %set nodes in x (nodes in y are identical)
    ds = mean(diff(x));    
    dt = ds^2/(alpha*4); %set timestep so as to ensure stability
    t(1,1,:) = 0:dt:tfinal;
%% Analytic Solution
    n = 20; %set number of retained eigenvalues
    Bi = h*L/k;
    Ta = zeros(N,N,length(t));
    Xn = zeros(N,N,length(t),n);
    theta = Ta;
    lg = [pi/(L*4)  pi/L:pi/L:(n-1)*pi/L]; %initialize guess for eigenvalues
    ln = zeros(1,n);
    for i = 1:length(lg)
        ln(i) = fzero(@(x) L*x.*tan(L*x)-Bi,lg(i)); %find eigenvalues 
    end
    xmap = repmat(x,[length(x) 1 length(t)]);
    tmap = repmat(t,[length(x) length(x) 1]);
    for i = 1:n
       Xn(:,:,:,i) = 2*sin(ln(i)*L)/(ln(i)*(L+sin(2*ln(i)*L)/(2*ln(i))))*exp(-(ln(i)^2)*tmap/alpha).*cos(ln(i)*xmap); 
    end
    X = sum(Xn,4);
    for i = 1:length(t)
        theta(:,:,i) = X(:,:,i).*X(:,:,i)';
    end
    Ta = Tamb+(To-Tamb)*theta;
%% Numeric Solution       
    Tn = zeros(N,N,length(t));
    Tn(:,:,1) = To; %Set initial temperature at all nodes
    
    f = waitbar(1/length(t),'Iterating...');
    for l = 2:length(t) %set timestepping loop
        Tn(N,1:N-1,l-1) = (Tn(N-1,1:N-1,l-1)+h/(k*ds)*Tamb)/(1+h/(k*ds)); %set top boundary (convection)
        Tn(1:N-1,N,l-1) = (Tn(1:N-1,N-1,l-1)+h/(k*ds)*Tamb)/(1+h/(k*ds)); %set right boundary (convection)
        Tn(N,N,l-1) = mean([Tn(N-1,N,l-1) Tn(N,N-1,l-1)]);
        Tn(1,:,l-1) = Tn(2,:,l-1); %set bottom boundary (dT/dy = 0)
        Tn(:,1,l-1) = Tn(:,2,l-1); %set left boundary (dT/dx = 0)
        d2Tdx2 = NumericalDerivative(2,1,ds,Tn(:,:,l-1));
        d2Tdy2 = NumericalDerivative(2,1,ds,Tn(:,:,l-1)')';
        Tn(2:N-1,2:N-1,l) = Tn(2:N-1,2:N-1,l-1) + alpha*dt*(d2Tdx2(2:N-1,2:N-1) + d2Tdy2(2:N-1,2:N-1)); %perform time step using forward Euler
        waitbar(l/length(t),f,['Time = ',num2str(t(l))]);
    end
    close(f)
end