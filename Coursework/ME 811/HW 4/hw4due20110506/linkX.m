function [AO,AE,AW,AN,AS,S] = linkX(U,V,P,UR,UL,UT,UB,VT,X,Y,rho,mu)

%this function calculates the link coefficients for the X-momentum equation
%given a guess for the velocity and pressure fields, and the mesh. Note
%that this is the specifically for a 2D staggered grid with the U-cell
%centered on the western boundary of the pressure-cell and the V-cell
%centered on the southern boundary. The code returns the link coefficients
%and source term for the X momentum equation.


%determining grid information
dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);
[N,M] = size(X);


%intializing matrices
AO = zeros(N,M);
AE = AO;
AW = AO;
AN = AO;
AS = AO;
S = AO;

Re = zeros(N,M);
Rw = Re;
Rn = Re;
Rs = Re;
MuN = Re;
MuS = Re;


%interpolating the values
%interior
for i = 2:M-1
    for j = 2:N-1
        Re(j,i) = rho(j,i)*(U(j,i+1) + U(j,i))/2;
        Rw(j,i) = rho(j,i-1)*(U(j,i) + U(j,i-1))/2;
        Rn(j,i) = (rho(j,i) + rho(j,i-1) + rho(j+1,i) + rho(j+1,i-1))/4*(V(j+1,i) + V(j+1,i-1))/2;
        Rs(j,i) = (rho(j,i) + rho(j,i-1) + rho(j-1,i) + rho(j-1,i-1))/4*(V(j,i) + V(j,i-1))/2;
        MuN(j,i) = (mu(j,i) + mu(j,i-1) + mu(j+1,i) + mu(j+1,i-1))/4;
        MuS(j,i) = (mu(j,i) + mu(j,i-1) + mu(j-1,i) + mu(j-1,i-1))/4;
    end
end

%right boundary
for j = 2:N-1
    Re(j,M) = rho(j,M)*(U(j,M)+UR(j))/2;
    Rw(j,M) = rho(j,M-1)*(U(j,M) + U(j,M-1))/2;
    Rn(j,M) = (rho(j,M) + rho(j,M-1) + rho(j+1,M) + rho(j+1,M-1))/4*(V(j+1,M) + V(j+1,M-1))/2;
    Rs(j,M) = (rho(j,M) + rho(j,M-1) + rho(j-1,M) + rho(j-1,M-1))/4*(V(j,M) + V(j,M-1))/2;
    MuN(j,M) = (mu(j,M) + mu(j,M-1) + mu(j+1,M) + mu(j+1,M-1))/4;
    MuS(j,M) = (mu(j,M) + mu(j,M-1) + mu(j-1,M) + mu(j-1,M-1))/4;
end

%top/bottom boundaries
for i = 2:M-1
    %top
    Re(N,i) = rho(N,i)*(U(N,i+1) + U(N,i))/2;
    Rw(N,i) = rho(N,i-1)*(U(N,i) + U(N,i-1))/2;
    Rn(N,i) = (rho(N,i) + rho(N,i-1))/2*(VT(i) + VT(i-1))/2;
    Rs(N,i) = (rho(N,i) + rho(N,i-1) + rho(N-1,i) + rho(N-1,i-1))/4*(V(N,i) + V(N,i-1))/2;
    MuN(N,i) = (mu(N,i) + mu(N,i-1))/2;
    MuS(N,i) = (mu(N,i) + mu(N,i-1) + mu(N-1,i) + mu(N-1,i-1))/4;
    
    %bottom
    Re(1,i) = rho(1,i)*(U(1,i+1) + U(1,i))/2;
    Rw(1,i) = rho(1,i-1)*(U(1,i) + U(1,i-1))/2;
    Rn(1,i) = (rho(1,i) + rho(1,i-1) + rho(2,i) + rho(2,i-1))/4*(V(2,i) + V(2,i-1))/2;
    Rs(1,i) = (rho(1,i) + rho(1,i-1))/2*(V(1,i) + V(1,i-1))/2;
    MuN(1,i) = (mu(1,i) + mu(1,i-1) + mu(2,i) + mu(2,i-1))/4;
    MuS(1,i) = (mu(1,i) + mu(1,i-1))/2;
end


%calculating the link coefficients
%interior
for i = 2:M-1
    for j = 2:N-1
        AO(j,i) = ((abs(Re(j,i)) + Re(j,i))/2 + (abs(Rw(j,i)) - Rw(j,i))/2)*dy + ...
                  ((abs(Rn(j,i)) + Rn(j,i))/2 + (abs(Rs(j,i)) - Rs(j,i))/2)*dx + ...
                  (mu(j,i) + mu(j,i-1))*dy/dx + (MuN(j,i) + MuS(j,i))*dx/dy;
        AE(j,i) = -(abs(Re(j,i)) - Re(j,i))/2*dy - mu(j,i)*dy/dx;
        AW(j,i) = -(abs(Rw(j,i)) + Rw(j,i))/2*dy - mu(j,i-1)*dy/dx;
        AN(j,i) = -(abs(Rn(j,i)) - Rn(j,i))/2*dx - MuN(j,i)*dx/dy;
        AS(j,i) = -(abs(Rs(j,i)) + Rs(j,i))/2*dx - MuS(j,i)*dx/dy;
        S(j,i) = dy*(P(j,i-1) - P(j,i));
    end
end


%right boundary
i = M;
for j = 2:N-1
    AO(j,i) = ((abs(Re(j,i)) + Re(j,i))/2 + (abs(Rw(j,i)) - Rw(j,i))/2)*dy + ...
              ((abs(Rn(j,i)) + Rn(j,i))/2 + (abs(Rs(j,i)) - Rs(j,i))/2)*dx + ...
              (mu(j,i) + mu(j,i-1))*dy/dx + (MuN(j,i) + MuS(j,i))*dx/dy;
    AE(j,i) = -(abs(Re(j,i)) - Re(j,i))/2*dy - mu(j,i)*dy/dx;
    AW(j,i) = -(abs(Rw(j,i)) + Rw(j,i))/2*dy - mu(j,i-1)*dy/dx;
    AN(j,i) = -(abs(Rn(j,i)) - Rn(j,i))/2*dx - MuN(j,i)*dx/dy;
    AS(j,i) = -(abs(Rs(j,i)) + Rs(j,i))/2*dx - MuS(j,i)*dx/dy;
    S(j,i) = dy*(P(j,i-1) - P(j,i)) - AE(j,i)*UR(j);
    AE(j,i) = 0;
end

%left boundary/left corners
%note all coefficients and source term have been initialized to zero
AO(:,1) = ones(N,1);
S(:,1) = UL;

%top and bottom boundaries
for i = 2:M-1
    %top
    AO(N,i) = ((abs(Re(N,i)) + Re(N,i))/2 + (abs(Rw(N,i)) - Rw(N,i))/2)*dy + ...
              (abs(Rs(N,i)) - Rs(N,i))/2*dx + ...
              (mu(N,i) + mu(N,i-1))*dy/dx + (3*MuN(N,i) + MuS(N,i))*dx/dy;
    AE(N,i) = -(abs(Re(N,i)) - Re(N,i))/2*dy - mu(N,i)*dy/dx;
    AW(N,i) = -(abs(Rw(N,i)) + Rw(N,i))/2*dy - mu(N,i-1)*dy/dx;
    AN(N,i) = 0;
    AS(N,i) = -(abs(Rs(N,i)) + Rs(N,i))/2*dx - (MuS(N,i) + MuN(N,i)/3)*dx/dy;
    S(N,i) = dy*(P(N,i-1) - P(N,i)) + (8/3*MuN(N,i)*dx/dy - Rn(N,i)*dx)*UT(i);
    
    %bottom
    AO(1,i) = ((abs(Re(1,i)) + Re(1,i))/2 + (abs(Rw(1,i)) - Rw(1,i))/2)*dy + ...
              (abs(Rn(1,i)) + Rn(1,i))/2*dx + ...
              (mu(1,i) + mu(1,i-1))*dy/dx + (MuN(1,i) + 3*MuS(1,i))*dx/dy;
    AE(1,i) = -(abs(Re(1,i)) - Re(1,i))/2*dy - mu(1,i)*dy/dx;
    AW(1,i) = -(abs(Rw(1,i)) + Rw(1,i))/2*dy - mu(1,i-1)*dy/dx;
    AN(1,i) = -(abs(Rn(1,i)) - Rn(1,i))/2*dx - (MuN(1,i) + MuS(1,i)/3)*dx/dy;
    AS(1,i) = 0;
    S(1,i) = dy*(P(1,i-1) - P(1,i)) + (8/3*MuS(1,i)*dx/dy + Rs(1,i)*dx)*UB(i);
end


%corners
%top-right
Re(N,M) = rho(N,M)*(U(N,M) + UR(N))/2;
Rw(N,M) = rho(N,M-1)*(U(N,M) + U(N,M-1))/2;
Rn(N,M) = (rho(N,M) + rho(N,M-1))/2*(VT(M) + VT(M-1))/2;
Rs(N,M) = (rho(N,M) + rho(N,M-1) + rho(N-1,M) + rho(N-1,M-1))/4*(V(N,M) + V(N,M-1))/2;
MuN(N,M) = (mu(N,M) + mu(N,M-1))/2;
MuS(N,M) = (mu(N,M) + mu(N,M-1) + mu(N-1,M) + mu(N-1,M-1))/4;

AO(N,M) = ((abs(Re(N,M)) + Re(N,M))/2 + (abs(Rw(N,M)) - Rw(N,M))/2)*dy + ...
          (abs(Rs(N,M)) - Rs(N,M))/2*dx + ...
          (mu(N,M) + mu(N,M-1))*dy/dx + (3*MuN(N,M) + MuS(N,M))*dx/dy;
AE(N,M) = -(abs(Re(N,M)) - Re(N,M))/2*dy - mu(N,M)*dy/dx;
AW(N,M) = -(abs(Rw(N,M)) + Rw(N,M))/2*dy - mu(N,M-1)*dy/dx;
AN(N,M) = 0;
AS(N,M) = -(abs(Rs(N,M)) + Rs(N,M))/2*dx - (MuS(N,M) + MuN(N,M)/3)*dx/dy;
S(N,M) = dy*(P(N,M-1) - P(N,M)) + (8/3*MuN(N,M)*dx/dy - Rn(N,M)*dx)*UT(M) - AE(N,M)*UR(N);
AE(N,M) = 0;

%bottom-right
Re(1,M) = rho(1,M)*(U(1,M) + UR(1))/2;
Rw(1,M) = rho(1,M-1)*(U(1,M) + U(1,M-1))/2;
Rn(1,M) = (rho(1,M) + rho(1,M-1) + rho(2,M) + rho(2,M-1))/4*(V(2,M) + V(2,M-1))/2;
Rs(1,M) = (rho(1,M) + rho(1,M-1))/2*(V(1,M) + V(2,M-1))/2;
MuN(1,M) = (mu(1,M) + mu(1,M-1) + mu(2,M) + mu(2,M-1))/4;
MuS(1,M) = (mu(1,M) + mu(1,M-1))/2;

AO(1,M) = ((abs(Re(1,M)) + Re(1,M))/2 + (abs(Rw(1,M)) - Rw(1,M))/2)*dy + ...
          (abs(Rn(1,M)) + Rn(1,M))/2*dx + ...
          (mu(1,M) + mu(1,M-1))*dy/dx + (MuN(1,M) + 3*MuS(1,M))*dx/dy;
AE(1,M) = -(abs(Re(1,M)) - Re(1,M))/2*dy - mu(1,M)*dy/dx;
AW(1,M) = -(abs(Rw(1,M)) + Rw(1,M))/2*dy - mu(1,M-1)*dy/dx;
AN(1,M) = -(abs(Rn(1,M)) - Rn(1,M))/2*dx - (MuN(1,M) + MuS(1,M)/3)*dx/dy;
AS(1,M) = 0;
S(1,M) = dy*(P(1,M-1) - P(1,M)) + (8/3*MuS(1,M)*dx/dy + Rs(1,M)*dx)*UB(M) - AE(1,M)*UR(1);
AE(1,M) = 0;