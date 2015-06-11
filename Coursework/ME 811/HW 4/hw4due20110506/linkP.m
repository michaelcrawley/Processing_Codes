function [AO,AE,AW,AN,AS,S] = linkP(U,V,Ax,Ay,X,Y,rho)

%this function calculates the link coefficients for the pressure correction
%equation given the x and y grid matrices, the central cell coefficients
%for the x and y momentum equations, the predicted velocity field, and the
%density field
%Note that this is the specifically for a 2D staggered grid with the U-cell
%centered on the western boundary of the pressure-cell and the V-cell
%centered on the southern boundary. The code returns the link coefficients
%and source term for the pressure correction equation.


%determining the grid size and initializing the variables
[N,M] = size(X);
dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

AO = zeros(N,M);
AE = AO;
AW = AO;
AN = AO;
AS = AO;
S = AO;

Re = zeros(N+2,M+2);
Rw = Re;
Rn = Re;
Rs = Re;

tempR = Re;
tempR(2:N+1,2:M+1) = rho;
rho = tempR;


%interpolating the density field
for i = 2:M+1
	for j = 2:N+1
		Re(j,i) = (rho(j,i+1) + rho(j,i))/2;
		Rw(j,i) = (rho(j,i-1) + rho(j,i))/2;
		Rn(j,i) = (rho(j+1,i) + rho(j,i))/2;
		Rs(j,i) = (rho(j-1,i) + rho(j,i))/2;
	end
end
Re = Re(2:N+1,2:M+1)*2;
Re(1:N,1:M-1) = Re(1:N,1:M-1)/2;

Rw = Rw(2:N+1,2:M+1)*2;
Rw(1:N,2:M) = Rw(1:N,2:M)/2;

Rn = Rn(2:N+1,2:M+1)*2;
Rn(1:N-1,1:M) = Rn(1:N-1,1:M)/2;

Rs = Rs(2:N+1,2:M+1)*2;
Rs(2:N,1:M) = Rs(2:N,1:M)/2;


%calculating the non-central cell coefficients
%interior
for i = 2:M-1
	for j = 2:N-1
		AE(j,i) = -Re(j,i)*dy^2/Ax(j,i+1);
		AW(j,i) = -Rw(j,i)*dy^2/Ax(j,i);
		AN(j,i) = -Rn(j,i)*dx^2/Ay(j+1,i);
		AS(j,i) = -Rs(j,i)*dx^2/Ay(j,i);
		S(j,i) = dy*(Rw(j,i)*U(j,i) - Re(j,i)*U(j,i+1)) + dx*(Rs(j,i)*V(j,i) - Rn(j,i)*V(j+1,i));
	end
end

%boundaries
%east
for j = 2:N-1
	AE(j,M) = 0;
	AW(j,M) = -Rw(j,M)*dy^2/Ax(j,M);
	AN(j,M) = -Rn(j,M)*dx^2/Ay(j+1,M);
	AS(j,M) = -Rs(j,M)*dx^2/Ay(j,M);
	S(j,M) = dy*Rw(j,M)*U(j,M) + dx*(Rs(j,M)*V(j,M) - Rn(j,M)*V(j+1,M));
end

%west
for j = 2:N-1
	AE(j,1) = -Re(j,1)*dy^2/Ax(j,2);
	AW(j,1) = 0;
	AN(j,1) = -Rn(j,1)*dx^2/Ay(j+1,1);
	AS(j,1) = -Rs(j,1)*dx^2/Ay(j,1);
	S(j,1) = -dy*Re(j,1)*U(j,2) + dx*(Rs(j,1)*V(j,1) - Rn(j,1)*V(j+1,1));
end

%north
for i = 2:M-1
	AE(N,i) = -Re(N,i)*dy^2/Ax(N,i+1);
	AW(N,i) = -Rw(N,i)*dy^2/Ax(N,i);
	AN(N,i) = 0;
	AS(N,i) = -Rs(N,i)*dx^2/Ay(N,i);
	S(N,i) = dy*(Rw(N,i)*U(N,i) - Re(N,i)*U(N,i+1)) + dx*Rs(N,i)*V(N,i);
end

%south
for i = 2:M-1
	AE(1,i) = -Re(1,i)*dy^2/Ax(1,i+1);
	AW(1,i) = -Rw(1,i)*dy^2/Ax(1,i);
	AN(1,i) = -Rn(1,i)*dx^2/Ay(2,i);
	AS(1,i) = 0;
	S(1,i) = dy*(Rw(1,i)*U(1,i) - Re(1,i)*U(1,i+1)) - dx*Rn(1,i)*V(2,i);
end


%corners
%north-east
AE(N,M) = 0;
AW(N,M) = -Rw(N,M)*dy^2/Ax(N,M);
AN(N,M) = 0;
AS(N,M) = -Rs(N,M)*dx^2/Ay(N,M);
S(N,M) = dy*Rw(N,M)*U(N,M) + dx*Rs(N,M)*V(N,M);

%north-west
AE(N,1) = -Re(N,M)*dy^2/Ax(N,2);
AW(N,1) = 0;
AN(N,1) = 0;
AS(N,1) = -Rs(N,M)*dx^2/Ay(N,1);
S(N,1) = -dy*Re(N,M)*U(N,2) + dx*Rs(N,M)*V(N,1);

%south-east
AE(1,M) = 0;
AW(1,M) = -Rw(1,M)*dy^2/Ax(1,M);
AN(1,M) = -Rn(1,M)*dx^2/Ay(2,M);
AS(1,M) = 0;
S(1,M) = dy*Rw(1,M)*U(1,M) - dx*Rn(1,M)*V(2,M);

%south-west
AE(1,1) = -Re(1,1)*dy^2/Ax(1,2);
AW(1,1) = 0;
AN(1,1) = -Rn(1,1)*dx^2/Ay(2,1);
AS(1,1) = 0;
S(1,1) = -dy*Re(1,1)*U(1,2) - dx*Rn(1,1)*V(2,1);


%calculating the central cell coefficients
AO = -(AE + AW + AN + AS);