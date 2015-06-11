%Completed by Michael Crawley for AAE 694 HW#1

Prms1 = 0.2; 
Prms2 = 0.2518;
w1 = 125*2*pi;
w2 = 160*2*pi;
t = 0:0.0001:0.2;
plot(t,Prms1*cos(w1*t)+Prms2*cos(w2*t));xlabel('t (s)');ylabel('P_t (Pa)');title('Sound Pressure level of combined wave');