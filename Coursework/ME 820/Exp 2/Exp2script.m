clear;clc;

Te = 21+273;
c = sqrt(1.4*287*Te);

%%Section 2
tmp = xlsread('experiment2 data.xls','Sheet1','A7:G402');
h=figure;plot(tmp(:,1),tmp(:,2));xlabel('frequency (Hz)');ylabel('Transmission Loss');title('Quarter Wave Resonator Transmission loss');xlim([48 1600]);
saveas(h,'QWR transmission loss','fig');saveas(h,'QWR transmission loss','png');
close(h);
[~,fI] = max(tmp(:,2));
f = tmp(fI,1);
lest = c/(4*f);
fprintf('Section 2: The estimated length is %6.3f m, the actual length is %6.3f m \n',[lest 0.854]);

%%Section 3
Dp = 4.859/100; Ap = 0.25*pi*Dp^2;
h=figure;plot(tmp(:,1),tmp(:,4));xlabel('frequency (Hz)');ylabel('Transmission Loss');title('Expansion Transmission loss');
saveas(h,'Expansion transmission loss','fig');saveas(h,'Expansion transmission loss','png');
TLmax = input('Specify TLmax from plot: ');
f = input('Specify Frequency from plot: ');
m = fzero(@(m) TLmax-10*log10(1+0.25*(m-1/m)^2),2);
Dest = sqrt(4*m*Ap/pi);
lest = c/(4*f);
fprintf('Section 3: The estimated diameter is %6.6f m, the estimated length is %6.6f m \n',[Dest lest]);
close(h);
%%Section 4
Lshort = 6.67/100;
Llong = 25.72/100;
fshortest = c/2/Lshort;
flongest = c/2/Llong;

h=figure;plot(tmp(:,1),tmp(:,5));xlabel('frequency (Hz)');ylabel('Transmission Loss');title('Short Perforated Duct Transmission loss'); 
saveas(h,'Perf short transmission loss','fig');saveas(h,'Perf short transmission loss','png');
fshort = input('Specify Frequency from plot: ');
close(h);

h=figure;plot(tmp(:,1),tmp(:,6));xlabel('frequency (Hz)');ylabel('Transmission Loss');title('Long Perforated Duct Transmission loss'); 
saveas(h,'Perf long transmission loss','fig');saveas(h,'Perf long transmission loss','png');
flong = input('Specify Frequency from plot: ');
close(h);

fprintf('Section 4a: The estimated frequency is %6.0f Hz, the actual frequency is %6.0f Hz \n',[fshortest fshort]);
fprintf('Section 4b: The estimated frequency is %6.0f Hz, the actual frequency is %6.0f Hz \n',[flongest flong]);

%%Section 5
sigmaS = 0.037;
sigmaL = 0.02;
teff = (0.079+0.75*0.249)/100;
d1 = 5.08/100;
d2S = 7.62/100;
d2L = 10.15/100;

fHL = c/pi*sqrt(sigmaL*d1/teff/(d2L^2-d1^2));
fHS = c/pi*sqrt(sigmaS*d1/teff/(d2S^2-d1^2));

fprintf('Section 5a: The estimated Helmholtz frequency for the short resonator is %6.0f Hz\n',fHS); 
fprintf('Section 5b: The estimated Helmholtz frequency for the long resonator is %6.0f Hz\n',fHL); 

%%Section 6
h=figure;plot(tmp(:,1),tmp(:,7));xlabel('frequency (Hz)');ylabel('Transmission Loss');title('Perforated Absorbing Silencer Transmission loss'); 
saveas(h,'Perf abs transmission loss','fig');saveas(h,'Perf abs transmission loss','png');
close(h);

%%Section 1
l = 8.5/100;
D = 4.044/100;
A = pi*(D^2)/4;

tmp = xlsread('experiment2 data.xls','Sheet1','I33:L803');
[~,fI] = max(tmp(:,2));
f = tmp(fI,1);
Vest = (((2*pi*f/c)^2)*l/A)^-1;
Vact = 0.2442*pi*(.1532^2)/4;
fprintf('Section 1: The estimated Volume is %6.6f m^3, the actual volume is %6.6f m^3 \n',[Vest Vact]);

%%Section 7
h=figure;plot(tmp(:,1),tmp(:,2),'k',tmp(:,1),tmp(:,3),'k-',tmp(:,1),tmp(:,4),'k*');xlabel('frequency (Hz)');ylabel('Transmission Loss');title('Helmholtz Resonator'); legend('M = 0.0', 'M = 0.05', 'M = 0.10');
saveas(h,'Helmholtz transmission loss','fig');saveas(h,'Helmholtz transmission loss','png');