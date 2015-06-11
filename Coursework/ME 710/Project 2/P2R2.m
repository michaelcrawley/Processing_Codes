function [T Nuxn] = P2R2(M,N)
    %M: Nodes in radial direction
    %N: Nodes is axial direction
    tic;
    n = 1;
    
    ro = .01; %m
    Tw = 500;
    Tin = 300;
    k = 0.6; %W/m/k
    rho = 998; %kg/m3
    v = 10^-6; %m2/s
    cp = 4182; %J/kg/K
    md = 0.01; %kg/s
    
    alpha = k/(rho*cp);
    Pr = v/alpha;
    um = md/(rho*pi*ro^2);
    Re = um*2*ro/v;
    
    L = 0.2*Re*Pr*ro;
    
    r = (0:ro/(M-1):ro)';
    rs = r/ro;
    x = 0:L/(N-2):L;
    dx = mean(diff(x));
    x = [-dx x];
    xs = alpha*x/(2*um*ro^2);    
    dr = mean(diff(r));
    u = 2*um*(1-rs.^2);
    
    xa = [0.001 0.004 0.01 0.04 0.08 0.10 0.20];%analytic axial locations
    Nuxa = [12.80 8.03 6.00 4.17 3.77 3.71 3.66];%analytic Nux solution
    
    %Numeric solution
    Y = zeros(N*M,1);
    Y(1:M) = 1; %set inlet boundary condition
    X = sparse(N*M);
    T = zeros(M,N,2);
    h = zeros(2,length(x));
    Nuxn = zeros(2,length(x));
    for s = 0:1
        A = u/dx+2*alpha/dr^2+2*s*alpha/dx^2;
        B = -s*alpha/dx^2;
        C = -2*alpha/dr^2;
        D = -(u/dx+s*alpha/dx^2);
        E = -(alpha/dr^2+alpha./(2*dr*r));
        F = -(alpha/dr^2-alpha./(2*dr*r));
        Ds = -(u/dx+2*s*alpha/dx^2);
        for i = 1:N
           for j = 1:M
              if i == 1 %inlet nodes
                  X(M*(i-1)+j,M*(i-1)+j) = 1;
              elseif i == 2 && j >= 2 && j <= M-2 %left boundary nodes
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i)+j) = B;
                  X(M*(i-1)+j,M*(i-1)+j+1) = E(j);
                  X(M*(i-1)+j,M*(i-1)+j-1) = F(j);
                  X(M*(i-1)+j,j) = D(j);
              elseif i==2 && j ==1 %bottom corner node
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i)+j) = B;
                  X(M*(i-1)+j,M*(i-1)+j+1) = C;
                  X(M*(i-1)+j,j) = D(j);
              elseif i==2 && j == M-1 %top corner node
                  X(M*(i-1)+j,M*(i-1)+j) =  A(j);
                  X(M*(i-1)+j,M*(i)+j) = B;
                  X(M*(i-1)+j,M*(i-1)+j-1) = F(j);
                  X(M*(i-1)+j,j) = D(j);
                  X(M*(i-1)+j,M*(i-1)+j+1) = E(j);
              elseif i>=3 && i<=N-1 &&j==1 %bottom boundary (symmetry)
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i-2)+j) = D(j);
                  X(M*(i-1)+j,M*(i)+j) = B;
                  X(M*(i-1)+j,M*(i-1)+j+1) = C;
              elseif i == N && j == 1 %bottom right corner node
                  X(M*(i-1)+j,M*(i-1)+j) =  A(j);
                  X(M*(i-1)+j,M*(i-2)+j) = D(j);
                  X(M*(i-1)+j,M*(i-1)+j+1) = C;
              elseif i>=2 && i <=N && j==M %top boundary (wall) nodes
                  X(M*(i-1)+j,M*(i-1)+j) = 1;
              elseif i>=3 && i <=N-1 && j==M-1 %top boundary (flow) nodes
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i-2)+j) = D(j);
                  X(M*(i-1)+j,M*(i)+j) = B;
                  X(M*(i-1)+j,M*(i-1)+j-1) = F(j);
                  X(M*(i-1)+j,M*(i-1)+j+1) = E(j);
              elseif i==N && j == M-1 %top right corner node
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i-2)+j) = Ds(j);
                  X(M*(i-1)+j,M*(i-1)+j-1) = F(j);
                  X(M*(i-1)+j,M*(i-1)+j+1) = E(j);
              elseif i==N && j>=2 && j<= M-2 %right boundary nodes
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i-2)+j) = Ds(j);
                  X(M*(i-1)+j,M*(i-1)+j-1) = F(j);
                  X(M*(i-1)+j,M*(i-1)+j+1) = E(j);
              else %interior nodes
                  X(M*(i-1)+j,M*(i-1)+j) = A(j);
                  X(M*(i-1)+j,M*(i-2)+j) = D(j);
                  X(M*(i-1)+j,M*(i)+j) = B;
                  X(M*(i-1)+j,M*(i-1)+j+1) = E(j);
                  X(M*(i-1)+j,M*(i-1)+j-1) = F(j);
              end
           end           
        end
        theta = reshape(X\Y,M,N);
        T(:,:,s+1) = Tw+(Tin-Tw)*theta;
        Tm = zeros(1,length(x));
        for i = 1:length(x)
           Tm(i) = 2/(um*ro^2)*trapz(r,u.*T(:,i,s+1).*r);
        end
        dTdr = NumericalDerivative(1,n,dr,T(end-2*n:end,:,s+1)')';
        h(s+1,:) = -k*dTdr(end,:)./(Tm-Tw);
        Nuxn(s+1,:) = 2*h(s+1,:)*ro/k;
    end
    ctime = toc;
    if exist(num2str(date),'file') ~= 7
        mkdir(pwd,num2str(date));
    end
    cd(num2str(date));
    filename = [num2str(M),' x ',num2str(N)];
    f(1) = figure;
    subplot(2,1,1),plot(xa,Nuxa,'*',xs(2:end),Nuxn(1,2:end));legend('Analytic','Numeric (axial conduction neglected)');title('Analytic vs Numeric Solution');ylabel('Nu_x');xlabel('x^*');
    subplot(2,1,2),plot(xa,100*(Nuxa-interp1(xs(2:end),Nuxn(1,2:end),xa,'spline'))./Nuxa),title('Nu_x (analytic) - Nu_x (numeric)');xlabel('x^*');ylabel('% error');
    f(2) = figure;
    subplot(2,1,2),plot(xs(2:end),Nuxn(2,2:end)-Nuxn(1,2:end));title('Nu_x (w/ axial conduction) - Nu_x (w/o axial conduction)');xlabel('x^*');ylabel('error');
    subplot(2,1,1),plot(xs(2:end),Nuxn(1,2:end),xs(2:end),Nuxn(2,2:end),'*');legend('Axial Conduction Neglected','Axial Conduction Included');title('Effect of Axial Conduction on Numeric Solution');xlabel('x^*');ylabel('Nu_x');
    f(3) = figure;
    contourf(x(2:end),r,(T(:,2:end,1)-Tw)/(Tin-Tw));colorbar; colormap gray; title('\theta');xlabel('x');ylabel('r');
    saveas(f(1),['Nuxa vs Nuxn ',filename],'fig');
    saveas(f(2),['Axial conduction effects ',filename],'fig');
    saveas(f(3),['Nondimensional temperature ',filename],'fig');
    close(f);
    save([filename,'.mat']);
    cd ..
end