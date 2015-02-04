%ME 804, Midterm,Problem 2
%Written by Michael Crawley
clear;clc;

%Gas Constants
R = 8.314; % in kJ/Kmole K
Na = 6.022e23; %Avogadro constant
k = 1.38e-23;
hbar = 1.05457E-34;

O2.y = 0.1;
O2.M = 32;
O2.Cp = 0.919*O2.M; %in kJ/Kmole
O2.theta = 2239; %Characteristic vibrational temperature, in K
O2.d = 3.9e-10;
O2.mu = 16e-3/Na;
O2.delta = 7.5e-3; %anharmonicity factor
O2.m = O2.mu*2;
O2.a = 0.451e-10;
O2.omega = O2.theta*k/hbar;

Ar.y = 0.9;
Ar.M = 40;
Ar.Cp = 0.52*Ar.M; %in kJ/Kmole

Cp = O2.y*O2.Cp+Ar.y*Ar.Cp;
gamma = Cp/(Cp-R);
Mix_M = O2.y*O2.M+Ar.y*Ar.M;

%Isentropic Flow calcs
As = 1; %throat area
Me = 8; %exit Mach number
Ae = ((gamma+1)/2)^(-(gamma+1)/(2*gamma-2))*(1+0.5*(gamma-1)*Me*Me)^((gamma+1)/(2*gamma-2))/Me; %exit Area ratio
N = 100; %number of points
x = linspace(0,5,N);
dx = mean(diff(x));
A = linspace(As,Ae,N); %calc area ratio at each point
To =  3000; % stagnation temperature, in K
Po = 100*101325; %stagnation pressure, in Pa
rho_o = Po*Mix_M/R/To;

M = zeros(1,N); %calc Mach number at each point
for j = 1:N
    M(j) = fzero(@(Mx) A(j)-((gamma+1)/2)^(-(gamma+1)/(2*gamma-2))*(1+0.5*(gamma-1)*Mx.*Mx).^((gamma+1)/(2*gamma-2))./Mx,[1 10]);
end

T = To*(1+0.5*(gamma-1)*M.^2).^-1;
P = Po*(1+0.5*(gamma-1)*M.^2).^(-gamma/(gamma-1));
rho = rho_o*(1+0.5*(gamma-1)*M.^2).^(-1/(gamma-1));
n = Na*rho/Mix_M;
c = sqrt(gamma*R*1000*T/Mix_M);
U = M.*c;
flow_time = dx./U;
% Ueq = O2.y*n*k*O2.theta*(0.5+1./(exp(O2.theta./T)-1)); %equilibrium vibrational energy for O2

Eeq = O2.y*n*k*O2.theta./((exp(O2.theta./T)-1)); %equilibrium vibrational energy for O2
Etheta =  O2.y*n*k/(exp(1)-1); 
kvvn1 = @(v,d,m,theta,thetap,delta,T) v*4*pi*(d^2)*sqrt(k*T/3/m).*(thetap/theta).*((thetap./T).^(1/6)).*exp(-1.5*(theta./T).^(1/3)).*exp(theta*(1-2*delta*v)/2./T)*exp(delta*v);
thetap = @(a,omega,mu) 4*pi*pi*omega*omega*mu*a*a/k;

k10 = kvvn1(1,O2.d,O2.m,O2.theta,thetap(O2.a,O2.omega,O2.mu),O2.delta,T);
tau = 1./(O2.y*n.*k10.*(1-exp(-O2.theta./T)));
% tau = k*T./(P.*k10.*(1-exp(-O2.theta./T)));
Ev = Eeq+(Etheta-Eeq).*exp(-flow_time./tau);
Tv = O2.theta./log(1+O2.y*n*k*O2.theta./Ev);

% Ev = zeros(1,N);
% Ev(1) = Eeq(1);
% for j = 2:N
%     Ev(j) = flow_time(j-1)*(Eeq(j)-Ev(j-1))/tau(j)+Ev(j-1);
% end
% Tv = O2.theta./log(1+O2.y*n*k*O2.theta./Ev); 

%plot figures
h(1) = figure;
plot(x,M); title('Mach number versus distance along nozzle'); xlabel('x (m)'); ylabel('Mach number');

h(2) = figure;
plot(x,P/1000); title('Pressure versus distance along nozzle'); xlabel('x (m)'); ylabel('Pressure (kPa)');

h(3) = figure;
plot(x,T); title('Temperature versus distance along nozzle'); xlabel('x (m)'); ylabel('Temp (K)');

h(4) = figure;
plot(x,n); title('Number density versus distance along nozzle'); xlabel('x (m)'); ylabel('n');

h(5) = figure;
plot(x,U); title('Flow velocity versus distance along nozzle'); xlabel('x (m)'); ylabel('U (m/s)');

h(6) = figure;
plot(x,flow_time); title('Flow time versus distance along nozzle'); xlabel('x (m)'); ylabel('t (s)');

h(7) = figure;
plot(x,T,x,Tv); title('Equilibrium temperature versus distance along nozzle'); xlabel('x (m)'); ylabel('T (k)'); legend('T','T_v');

h(8) = figure;
plot(x,tau); title('Vibrational relaxation time versus distance along nozzle'); xlabel('x (m)'); ylabel('\tau (s)'); 

for j = 1:length(h)
   saveas(h(j),num2str(j),'fig'); 
end

close(h);