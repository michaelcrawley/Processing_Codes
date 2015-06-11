%ME 804, Midterm,Problem 1
%Written by Michael Crawley
clear;clc;

%Constants
hbar = 1.05457E-34;
k = 1.3806e-23;
Na = 6.022e23; %Avogadro constant

%functions
kvvn1 = @(v,d,m,theta,thetap,delta,T) v*4*pi*(d^2)*sqrt(k*T/3/m)*(thetap/theta)*((thetap/T).^(1/6))*exp(-1.5*(theta/T).^(1/3))*exp(theta*(1-2*delta*v)/2/T)*exp(delta*v);
thetap = @(a,omega,mu) 4*pi*pi*omega*omega*mu*a*a/k;

%Experimental data
ptau = @(A,u,X) 101325*exp(A*(X-0.015*u^0.25)-18.42); %note that this is not T^(-1/3), units have been converted to Pa-s

N2.A = 220;
N2.theta = 3395; %in K
N2.u = 14;
N2.X = 0.03:0.01:0.09;
N2.T = N2.X.^-3;
N2.d = 2.9e-10;
N2.mu =  14e-3/Na; %in kg
N2.delta = 6e-3; % anharmonicity factor
N2.dis = 942000/Na/k; %dissociation energy in Kelvin

O2.A = 129;
O2.theta = 2239;
O2.u = 16;
O2.X = 0.05:0.01:0.12;
O2.T = O2.X.^-3;
O2.d = 3.9e-10;
O2.mu = 16e-3/Na;
O2.delta = 7.5e-3; %anharmonicity factor
O2.dis = 494000/Na/k; %dissociation energy in Kelvin

%Calculate a for N2
N2.omega = N2.theta*k/hbar;
N2.ptauex = ptau(N2.A,N2.mu,N2.X); 
N2.dE = hbar*N2.omega;
N2.m = N2.mu*2;
N2.RHS = (N2.theta*N2.T.^(2/3)*sqrt(3*N2.m*k).*exp(-N2.dE/2/k./N2.T))./(4*pi*(N2.d^2)*N2.ptauex.*(1-exp(-N2.theta./N2.T))); %define RHS of eq. 1
N2.fitfun = @(x,a) (thetap(a,N2.omega,N2.mu).^(7/6)).*exp(-1.5*(thetap(a,N2.omega,N2.mu)./x).^(1/3));

% N2.a = GeneralFit(N2.T,N2.RHS,N2.fitfun);
N2.a = 0.605e-10;
N2.vd = floor((-1/2/N2.delta/N2.theta)*(-(1-N2.delta)*N2.theta-sqrt(((1-N2.delta)*N2.theta)^2-4*(N2.delta*N2.theta*N2.dis))));
N2.k1 = kvvn1(1,N2.d,N2.m,N2.theta,thetap(N2.a,N2.omega,N2.mu),N2.delta,5000); %k10 rate at T = 5000 K
N2.kvd = kvvn1(N2.vd,N2.d,N2.m,N2.theta,thetap(N2.a,N2.omega,N2.mu),N2.delta,5000); %k vd->(vd-1) rate at T = 5000 K

%Calculate a for O2
O2.omega = O2.theta*k/hbar;
O2.ptauex = ptau(O2.A,O2.mu,O2.X); 
O2.dE = hbar*O2.omega;
O2.m = O2.mu*2;
O2.RHS = (O2.theta*O2.T.^(2/3)*sqrt(3*O2.m*k).*exp(-O2.dE/2/k./O2.T))./(4*pi*(O2.d^2)*O2.ptauex.*(1-exp(-O2.theta./O2.T))); %define RHS of eq. 1
O2.fitfun = @(x,a) (thetap(a,O2.omega,O2.mu).^(7/6)).*exp(-1.5*(thetap(a,O2.omega,O2.mu)./x).^(1/3));

% O2.a = GeneralFit(O2.T,O2.RHS,O2.fitfun);
O2.a = 0.451e-10;
O2.vd = floor((-1/2/O2.delta/O2.theta)*(-(1-O2.delta)*O2.theta-sqrt(((1-O2.delta)*O2.theta)^2-4*(O2.delta*O2.theta*O2.dis))));
O2.k1 = kvvn1(1,O2.d,O2.m,O2.theta,thetap(O2.a,O2.omega,O2.mu),O2.delta,5000); %k10 rate at T = 5000 K
O2.kvd = kvvn1(O2.vd,O2.d,O2.m,O2.theta,thetap(O2.a,O2.omega,O2.mu),O2.delta,5000); %k vd->(vd-1) rate at T = 5000 K

%Part D
T = 3140.7; %in K
p = 76000; %in Pa
N2.k1_N2 = kvvn1(1,N2.d,N2.m,N2.theta,thetap(N2.a,N2.omega,N2.mu),N2.delta,T); % N2-N2 reaction rate
N2.k1_O2 = kvvn1(1,N2.d,N2.m,N2.theta,thetap(N2.a,N2.omega,(28*32)/(28+32)/1000/Na),N2.delta,T); % N2-O2 reaction rate
N2.k1t = 0.79*N2.k1_N2+0.21*N2.k1_O2;
N2.tau = k*T/(p*N2.k1t*(1-exp(-N2.theta/T)));

O2.k1_O2 = kvvn1(1,O2.d,O2.m,O2.theta,thetap(O2.a,O2.omega,O2.mu),O2.delta,T); % O2-O2 reaction rate
O2.k1_N2 = kvvn1(1,O2.d,O2.m,O2.theta,thetap(O2.a,O2.omega,(28*32)/(28+32)/1000/Na),O2.delta,T); % O2-N2 reaction rate
O2.k1t = 0.79*O2.k1_N2+0.21*O2.k1_O2;
O2.tau = k*T/(p*O2.k1t*(1-exp(-O2.theta/T)));