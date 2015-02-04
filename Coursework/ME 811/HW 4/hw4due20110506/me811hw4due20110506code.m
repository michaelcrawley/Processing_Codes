%ME811 HW4 Due 5-6-2011 Code   Nathan Webb
close all
clear
clc


%constants
M = 5; %number of cells in the x direction
N = 5; %number of cells in the y direction

rho = 1000; %kg/m^3 density
mu = 0.001; %kg/(m*s) viscosity
Uo = 0.01; %m/s lid velocity
L = 0.01; %m cavity dimension (square cavity)

Wuv = 0.5; %linear relaxation for the velocity fields
Wp = 0.2; %linear relaxation for the pressure field
alpha = 0.2; %inertial relaxation for the velocity fields

UVsweeps = 2; %number of sweeps to perform for the X/Y-momentum equations
Psweeps = 20; %number of sweeps to perform for the pressure correction equation
maxIT = 10; %maximum number of outer iterations

UVtol = 10^-10; %residual tolerance for the momentum equation solution
Ptol = 10^-10; %residual tolerance for the pressure correction equation

fname = 'HighRe.mat';



%INSERT GRID LOOP HERE


%calculating the grid
[Xp,Yp,Xu,Yu,Xv,Yv] = staggered_grid(L,L,M,N);
Rho = rho*ones(size(Xp));
Mu = mu*ones(size(Xp));

dx = Xp(1,2) - Xp(1,1);
dy = Yp(2,1) - Yp(1,1);


%INSERT REYNOLDS LOOP HERE


%initially guessing the velocity/pressure fields
U = zeros(size(Xu));
V = zeros(size(Xv));
P = zeros(size(Xp));

UT = Uo*ones(1,M);
UB = 0*UT;
UL = zeros(N,1);
UR = UL;

VT = zeros(1,M);
VB = VT;
VL = zeros(N,1);
VR = VL;

U(:,1) = UL;
V(1,:) = VB;


%calculating the initial residual
[AxO,AxE,AxW,AxN,AxS,Sx] = linkX(U,V,P,UR,UL,UT,UB,VT,Xu,Yu,Rho,Mu);
[AyO,AyE,AyW,AyN,AyS,Sy] = linkY(U,V,P,UR,VR,VL,VT,VB,Xv,Yv,Rho,Mu);
[Ru,Rv,Riu,Riv] = residXY(U,AxO,AxE,AxW,AxN,AxS,Sx,V,AyO,AyE,AyW,AyN,AyS,Sy);

[ApO,ApE,ApW,ApN,ApS,Sp] = linkP(U,V,AxO,AyO,Xp,Yp,Rho);
Rp = residP(P,ApO,ApE,ApW,ApN,ApS,Sp);
    

%outer iteration loop
n = 1;
while ((Ru(n) >= UVtol) || (Rv(n) >= UVtol) || (Rp(n) >= Ptol)) && (n <= maxIT)

    %solving the X/Y-momentum equations for u-hat and v-hat
    uP = solveEQ((1 + alpha)*AxO,AxE,AxW,AxN,AxS,Riu,zeros(size(Xu)),UVsweeps);
    vP = solveEQ((1 + alpha)*AyO,AyE,AyW,AyN,AyS,Riv,zeros(size(Xv)),UVsweeps);
	
	
	%calculating u-hat and v-hat
	uH = U + uP;
	vH = V + vP;


    %calculating the pressure link coefficients
    [ApO,ApE,ApW,ApN,ApS,Sp] = linkP(uH,vH,AxO,AyO,Xp,Yp,Rho);


    %solving the pressure correction equation
    Pp = solveEQ(ApO,ApE,ApW,ApN,ApS,Sp,zeros(size(Xp)),Psweeps);


    %updating the velocity and pressure fields
	uPP = dy*([zeros(N,1),Pp(:,1:M-1)] - Pp)./AxO;
	vPP = dx*([zeros(1,M);Pp(1:N-1,:)] - Pp)./AyO;
    
    uPP(:,1) = 0*UL;
    vPP(1,:) = 0*VB;
	
    U = uH + Wuv*uPP;
    V = vH + Wuv*vPP;
	
	P = P + Wp*Pp;
    
    
    %calculating the residuals
    [AxO,AxE,AxW,AxN,AxS,Sx] = linkX(U,V,P,UR,UL,UT,UB,VT,Xu,Yu,Rho,Mu);
    [AyO,AyE,AyW,AyN,AyS,Sy] = linkY(U,V,P,UR,VR,VL,VT,VB,Xv,Yv,Rho,Mu);
    [Ru(n+1),Rv(n+1),Riu,Riv] = residXY(U,AxO,AxE,AxW,AxN,AxS,Sx,V,AyO,AyE,AyW,AyN,AyS,Sy);

    [ApO,ApE,ApW,ApN,ApS,Sp] = linkP(U,V,AxO,AyO,Xp,Yp,Rho);
    Rp(n+1) = residP(Pp,ApO,ApE,ApW,ApN,ApS,Sp);
    

    n = n + 1;
	if rem(n,100) == 0
		fprintf('%d\n',n);
	end
    
end


%post processing residuals
Ru = Ru(2:end);
Rv = Rv(2:end);
Rp = Rp(2:end);

IT = [1:n-1];


%%%%%%%%%%%%%%%%%%%%%%% TEMPORARY %%%%%%%%%%%%%%%%%%%%%%%%
save(fname,'Xp','Yp','U','V','P','Ru','Rv','Rp');

[U,V] = stag_interp(U,V,UR,VT);
Ur = U(1:N/20:N,1:M/20:M);
Vr = V(1:N/20:N,1:M/20:M);
Xr = Xp(1:N/20:N,1:M/20:M);
Yr = Yp(1:N/20:N,1:M/20:M);

figure
quiver(Xr,Yr,Ur,Vr);
xlim([0,0.01]);
ylim([0,0.01]);

figure
contourf(Xp,Yp,P,15);
xlim([0,0.01]);
ylim([0,0.01]);
colorbar

figure
semilogy(IT,Ru,'-k',IT,Rv,'--k',IT,Rp,':k');
legend('U','V','P');
set(legend,'Location','NorthEast');