x = 1:20;
t = 2:7;
s = 1:10;
c = .1:0.1:2;

[X,T,S,C] = ndgrid(x,t,s,c);

fun1 = X.*cosd(X.*T);%this will represent the physical signal
fun2 = log(X + T.*S.*C); %this will represent the daughter wavelets

%Method 1: Multiply fun1 and fun2 first, then integrate over all dimensions
tmp = fun1.*fun2;
tmp = trapz(c,tmp,4);
tmp = trapz(s,tmp,3);
tmp = trapz(t,tmp,2);
Method1 = trapz(x,tmp);

%Method 2: Integrate fun2 over S,C then multiply with fun1, then integrate
%over remaining dimensions
tmp = trapz(c,fun2,4);
tmp = trapz(s,tmp,3);
tmp = fun1(:,:,1,1).*tmp;
tmp = trapz(t,tmp,2);
Method2 = trapz(x,tmp);