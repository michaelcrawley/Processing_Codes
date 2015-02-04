clear;
D = 1; %Code currently only works for first order derivatives
N = 3;
M = 3;
h = 1/2;
A = zeros(N+M+1,N+M+1); %initialize left matrix for solution (Ax=F)
F = zeros(N+M+1,1); %initialize right matrix for solution
F(D+1) = factorial(D); %Set derivative dependence 
syms b;

%Set N:N-2 rows in left matrix based off of taylor series expansions
for i=1:N+M
    A(i,:) = (-N:M).^(i-1);
end

%Set second to last row based off of left-right taylor series property
if N == M
    A(N+M,:) = [(-1)^(D+1)*ones(1,N) 0 ones(1,M)]; %This is only correct for central differencing
end

%calculate d/dan coefficients
da = [1;A(1:N+M,2:N+M+1)\-A(1:N+M,1)];

%Create integration string
str = ['(',num2str(da(1)),'*i*exp(i*',num2str(-N),'*b)'];
for i = -N+1:M
    str = strcat(str,['+(',num2str(da(i+N+1)),')*i*exp(i*',num2str(i),'*b)']);
end
str = [str,')'];

%Integrate for final linear equation
Y = double(int(['-b*',str],-pi/2,pi/2)); %Why does this need to be -b instead of b?
F(end) = Y;
for i = -N:M
    A(N+M+1,i+N+1)= double(int(['i*exp(i*',num2str(i),'*b)*',str],-pi/2,pi/2));
end

%Calculate coefficients from Ax=F
coefs = ((A\F)')/(h^(N+M-2));