function [Tn Ta] = P1R3(L, N, To, Tamb,tfinal)
%% Variable initialization
    cp = 0.485; 
    rho = 7.82E6;
    h= 100;
    k = 20;
    
    alpha = k/(rho*cp);
    x = 0:L/(N-1):L; %set nodes in x (nodes in y are identical)
    ds = mean(diff(x));    
    dt = ds^2/(alpha*10); %set timestep so as to ensure stability    
    t(1,1,:) = 0:dt:tfinal;
    
    Ta = zeros(N,N,length(t));
    n = 30; %set number of retained eigenvalues
    Xn = zeros(N,N,length(t),n);
    xmap = repmat(x,[length(x) 1 length(t)]);
    tmap = repmat(t,[length(x) length(x) 1]);
    
    Tn = zeros(N,N,length(t));
%% Analytic Solution    
    Bi = h*L/k;   
    theta = Ta;
    lg = [pi/(L*4)  pi/L:pi/L:(n-1)*pi/L]; %initialize guess for eigenvalues
    lambda = zeros(1,n);
    for i = 1:length(lg)
        lambda(i) = fzero(@(x) L*x.*sin(L*x)-Bi*cos(L*x),lg(i)); %find eigenvalues 
    end
    for i = 1:n
       Xn(:,:,:,i) = 2*sin(lambda(i)*L)/(lambda(i)*(L+sin(2*lambda(i)*L)/(2*lambda(i))))*exp(-(lambda(i)^2)*tmap*alpha).*cos(lambda(i)*xmap); 
    end
    X = sum(Xn,4);
    for i = 1:length(t)
        theta(:,:,i) = X(:,:,i).*X(:,:,i)';
    end    
    Ta = Tamb+(To-Tamb)*theta;
    Ta = [Ta(end:-1:2,end:-1:2,:) Ta(end:-1:2,:,:); Ta(:,end:-1:2,:) Ta(:,:,:)];
    if exist(num2str(date),'file') ~= 7
        mkdir(pwd,num2str(date));
    end
    cd(num2str(date));
    filename = [num2str(clock),'.mat'];
    save(filename, 'Ta','t','alpha','h','lambda', 'k','cp','rho','L','To','Tamb');
    clear X Xn theta Ta;
%% Numeric Solution 
    Tn(:,:,1) = To; %Set initial temperature at all nodes
    
    f = waitbar(1/length(t),'Iterating...');
    for l = 2:length(t) %set timestepping loop
        d2Tdx2 = NumericalDerivative(2,1,ds,Tn(:,:,l-1));
        d2Tdy2 = NumericalDerivative(2,1,ds,Tn(:,:,l-1)')';
        
        d2Tdy2(N,1:N-1) = 2*Tn(N-1,1:N-1,l-1)/(ds^2)-2*Tn(N,1:N-1,l-1)/(ds^2)+2*(-h/k*(Tn(N,1:N-1,l-1)-Tamb)/ds); %set right convection boundary condition
        d2Tdx2(1:N-1,N) = 2*Tn(1:N-1,N-1,l-1)/(ds^2)-2*Tn(1:N-1,N,l-1)/(ds^2)+2*(-h/k*(Tn(1:N-1,N,l-1)-Tamb)/ds); %set top convection boundary condition
        d2Tdy2(1,1:N) = 2*(Tn(2,1:N,l-1)-Tn(1,1:N,l-1))/(ds^2); %set left boundary condition (dT/dx = 0)
        d2Tdx2(1:N,1) = 2*(Tn(1:N,2,l-1)-Tn(1:N,1,l-1))/(ds^2); %set bottom boundary condition (dT/dy = 0);
        d2Tdy2(N,N) = 2*(Tn(N,N-1,l-1)-Tn(N,N,l-1))/(ds^2)+2*(-h/k*(Tn(N,N,l-1)-Tamb)/ds); %set corner boundary condition
        d2Tdx2(N,N) = 2*(Tn(N-1,N,l-1)-Tn(N,N,l-1))/(ds^2)+2*(-h/k*(Tn(N,N,l-1)-Tamb)/ds); %set corner boundary condition
        
        Tn(:,:,l) = Tn(:,:,l-1) + alpha*dt*(d2Tdx2 + d2Tdy2); %perform time step using forward Euler
        waitbar(l/length(t),f,['Time = ',num2str(t(l))]);
    end
    Tn = [Tn(end:-1:2,end:-1:2,:) Tn(end:-1:2,:,:); Tn(:,end:-1:2,:) Tn(:,:,:)];
    x = [-x(end:-1:2) x];
    close(f)
    save(filename, 'Tn', 'x','-append');
    cd ..
end