%Written by Michael Crawley for ME 818, HW #4, Problem #2
clear;clc;

fun = @(x) 1./(1+x.^2);
funhatanalytic = @(xi) pi*exp(-abs(xi));

N = 2^14;
L = 2^11;
omega = pi*N/2/L;

dx = 2*L/(N);
x = -L:dx:L-dx;
dxi = 2*pi/N/dx;
xi = -omega:dxi:omega-dxi;

f = fun(x);

k = 0:(N-1);
fhatnum = f*(dx*((-1).^k)'*((-1).^k).*exp(-2*pi*1i*(k'*k)/N));

% plot(xi,abs(fhatnum),xi,funhatanalytic(xi)); xlim([-10 10]);

I = N/2+1:N/2+41;

fid=fopen('H4P4.txt','w');
fprintf(fid,['ME 818 Homework 4, Problem 2\nCompleted by Michael Crawley on ',date,'\n']);
fprintf(fid,'N = %d, L = %d\n',[N; L]);
fprintf(fid,'\nxi\t\tNumeric Soln\tAnalytic Soln\tPercent Error\n');
fprintf(fid,'%1.3f\t\t %3.4f\t\t %3.4f\t\t %3.4f\n',[xi(I); fhatnum(I); funhatanalytic(xi(I));100*(fhatnum(I)-funhatanalytic(xi(I)))./funhatanalytic(xi(I)) ]);
fclose(fid);