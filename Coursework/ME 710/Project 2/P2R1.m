function [ftheta A F] = P2R1(N,M)
    %N: Nodes in radial direction
    %M: Nodes is axial direction

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
    
    r = (0:ro/(N-1):ro)';
    rs = r/ro;
    x = 0:L/(M-1):L;
    xs = alpha*x/(2*um*ro^2);
    dx = mean(diff(xs));
    dr = mean(diff(rs));
    u = 2*um*(1-rs.^2);
    
    xa = [0.001 0.004 0.01 0.04 0.08 0.10 0.20];%analytic axial locations
    Nuxa = [12.80 8.03 6.00 4.17 3.77 3.71 3.66];%analytic Nux solution
    
    %Numeric solution
    F = zeros(N*M,1);
    F(1:M) = 1; %set inlet boundary condition
    A = zeros(N*M);
    for i = 1:N
       for j = 1:M
          if i == 1 %inlet nodes
              A(M*(i-1)+j,M*(i-1)+j) = 1;
          elseif i == 2 && j >= 2 && j <= M-2 %left boundary nodes
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i)+j) = -alpha/dx^2;
              A(M*(i-1)+j,M*(i-1)+j+1) = -(alpha/dr^2+alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,M*(i-1)+j-1) = -(alpha/dr^2-alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,j) = -(u(j)/dx+alpha/dx^2);
          elseif i==2 && j ==1 %bottom corner node
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i)+j) = -alpha/dx^2;
              A(M*(i-1)+j,M*(i-1)+j+1) = -2*alpha/dr^2;
              A(M*(i-1)+j,j) = -(u(j)/dx+alpha/dx^2);
          elseif i==2 && j == M-1 %top corner node
              A(M*(i-1)+j,M*(i-1)+j) =  u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i)+j) = -alpha/dx^2;
              A(M*(i-1)+j,M*(i-1)+j-1) = -(alpha/dr^2-alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,j) = -(u(j)/dx+alpha/dx^2);
              A(M*(i-1)+j,M*(i-1)+j+1) = -(alpha/dr^2+alpha/(2*rs(j)*dr));
          elseif i>=3 && i<=N-1 &&j==1 %bottom boundary (symmetry)
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i-2)+j) = -(u(j)/dx+alpha/dx^2);
              A(M*(i-1)+j,M*(i)+j) = -alpha/dx^2;
              A(M*(i-1)+j,M*(i-1)+j+1) = -2*alpha/dr^2;
          elseif i == N && j == 1 %bottom right corner node
              A(M*(i-1)+j,M*(i-1)+j) =  u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i-2)+j) = -(u(j)/dx+2*alpha/dx^2);
              A(M*(i-1)+j,M*(i-1)+j+1) = -2*alpha/dr^2;
          elseif i>=2 && i <=N && j==M %top boundary (wall) nodes
              A(M*(i-1)+j,M*(i-1)+j) = 1;
          elseif i>=3 && i <=N-1 && j==M-1 %top boundary (flow) nodes
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i-2)+j) = -(u(j)/dx+alpha/dx^2);
              A(M*(i-1)+j,M*(i)+j) = -alpha/dx^2;
              A(M*(i-1)+j,M*(i-1)+j-1) = -(alpha/dr^2-alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,M*(i-1)+j+1) = -(alpha/dr^2+alpha/(2*rs(j)*dr));
          elseif i==N && j == M-1 %top right corner node
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i-2)+j) = -(u(j)/dx+2*alpha/dx^2);
              A(M*(i-1)+j,M*(i-1)+j-1) = -(alpha/dr^2-alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,M*(i-1)+j+1) = -(alpha/dr^2+alpha/(2*rs(j)*dr));
          elseif i==N && j>=2 && j<= M-2 %right boundary nodes
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i-2)+j) = -(u(j)/dx+2*alpha/dx^2);
              A(M*(i-1)+j,M*(i-1)+j-1) = -(alpha/dr^2-alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,M*(i-1)+j+1) = -(alpha/dr^2+alpha/(2*rs(j)*dr));
          else %interior nodes
              A(M*(i-1)+j,M*(i-1)+j) = u(j)/dx+2*alpha/dr^2+2*alpha/dx^2;
              A(M*(i-1)+j,M*(i-2)+j) = -(u(j)/dx+alpha/dx^2);
              A(M*(i-1)+j,M*(i)+j) = -alpha/dx^2;
              A(M*(i-1)+j,M*(i-1)+j+1) = -(alpha/dr^2+alpha/(2*rs(j)*dr));
              A(M*(i-1)+j,M*(i-1)+j-1) = -(alpha/dr^2-alpha/(2*rs(j)*dr));
          end
       end
    end
    
    theta = F\A;
    ftheta = zeros(N,M);
    for i = 1:N
       for j = 1:M
          ftheta(j,i) = theta(M*(i-1)+j); 
       end
    end
end