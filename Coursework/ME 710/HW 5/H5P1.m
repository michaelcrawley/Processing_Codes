%Problem 1, Homework 5, ME 710
T = 2000;
lambda = [0.28 0.36 0.48 1.25 2.6];
e = [0.43 0.47 0.47 0.32 0.18];

c = 2.998E8;
h = 6.626E-34;
k = 1.381E-23;
Ebl = 2*pi^2*h*c^2./(lambda.^5.*(exp(h*c./lambda/k/T)-1)); %this is not correct!

et = trapz(lambda,e.*Ebl)/trapz(lambda,Ebl);