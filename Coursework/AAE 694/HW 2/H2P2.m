%AAE 694, Homework 2, Problem 2
%Completed by Michael Crawley
clear;clc;
x{1} = 0:0.01:2;
x{2} = 0:0.2:2;

f = @(x)  x.^2-5*sin(pi*x);
dfdx = @(x) 2*x-5*pi*cos(pi*x);
d2fdx2 = @(x) 2+5*pi*pi*sin(pi*x);

method = {'-center', '-downstream' ,'-upstream'};
Horder = [1 2 2];

for i=1:3
    for j=1:2
        figure;
        subplot(4,1,1),plot(x{j},NumericalDerivative(1,Horder(i),mean(diff(x{j})),f(x{j}),method{i}));xlabel('x');ylabel('df/dx');%Plot first derivative
        title(['Derivative: 1st, Method: ', method{i}(2:end), ', Step Size: ',num2str(mean(diff(x{j})))]);
        subplot(4,1,2),plot(x{j},NumericalDerivative(1,Horder(i),mean(diff(x{j})),f(x{j}),method{i})-dfdx(x{j}));xlabel('x');ylabel('Numerical-Analytical');%Plot first derivative error
        title(['Error in Derivative: 1st, Method: ', method{i}(2:end), ', Step Size: ',num2str(mean(diff(x{j})))]);   
        subplot(4,1,3),plot(x{j},NumericalDerivative(2,Horder(i),mean(diff(x{j})),f(x{j}),method{i}));xlabel('x');ylabel('d2f/dx2');%Plot second derivative
        title(['Derivative: 2nd, Method: ', method{i}(2:end), ', Step Size: ',num2str(mean(diff(x{j})))]);
        subplot(4,1,4),plot(x{j},NumericalDerivative(2,Horder(i),mean(diff(x{j})),f(x{j}),method{i})-d2fdx2(x{j}));xlabel('x');ylabel('Numerical-Analytical');%Plot first derivative error
        title(['Error in Derivative: 2nd, Method: ', method{i}(2:end), ', Step Size: ',num2str(mean(diff(x{j})))]);
        figure;
        subplot(2,1,1),plot(x{j},NumericalDerivative(1,Horder(i),mean(diff(x{j})),f(x{j}),method{i}),x{j},dfdx(x{j}));
        subplot(2,1,2),plot(x{j},NumericalDerivative(2,Horder(i),mean(diff(x{j})),f(x{j}),method{i}),x{j},d2fdx2(x{j}));
    end
end