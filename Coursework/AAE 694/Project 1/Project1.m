% Project 1
% Completed by Michael Crawley for AAE 694
clear all;

tic;
%3rd order optimized time discretization coefficients (b0, b1, b2, b3)
b = [2.3025580883830 -2.4910075998482 1.5743409331815 -0.38589142217162];
dx = 1.0/2;
dt = 0.1/2;

%Case 1
x1 = -20:dx:450;
t1 = [0 0 0 0:dt:400];

%Initialize solutions
u2CDS1 = zeros(length(t1),length(x1)); %2nd Order Central Difference
u2CDS1(1,:) = 0.5*exp(-log(2)*(x1/3).^2);
u2CDS1(2,:) = u2CDS1(1,:);
u2CDS1(3,:) = u2CDS1(1,:);
u2CDS1(4,:) = u2CDS1(1,:);
u4CDS1 = u2CDS1; %4th Order Central Difference
u6CDS1 = u2CDS1; %6th Order Central Difference
u4DRP1 = u2CDS1; %4th Order DRP
ua1 = 0.5*exp(-log(2)*((x1-400)/3).^2); %analytic solution

for j = 4:length(t1)-1
    u2CDS1(j+1,:) = u2CDS1(j,:)-dt*b*NumericalDerivative(1,2,dx,u2CDS1(j:-1:j-3,:)); 
    u2CDS1(j+1,1) = 0;   %suppress left boundary fluctuation
    
    u4CDS1(j+1,:) = u4CDS1(j,:)-dt*b*NumericalDerivative(1,4,dx,u4CDS1(j:-1:j-3,:));
    u4CDS1(j+1,1:2) = 0;   %suppress left boundary fluctuation
    
    u6CDS1(j+1,:) = u6CDS1(j,:)-dt*b*NumericalDerivative(1,6,dx,u6CDS1(j:-1:j-3,:));
    u6CDS1(j+1,1:3) = 0;   %suppress left boundary fluctuation
    
    u4DRP1(j+1,:) = u4DRP1(j,:)-dt*b*FakeNumericalDerivative(dx,u4DRP1(j:-1:j-3,:));
    u4DRP1(j+1,1:3) = 0;   %suppress left boundary fluctuation    
end

%calc rms error
u2CDS1error = sqrt(sum((ua1-u2CDS1(end,:)).^2));
u4CDS1error = sqrt(sum((ua1-u4CDS1(end,:)).^2));
u6CDS1error = sqrt(sum((ua1-u6CDS1(end,:)).^2));
u4DRP1error = sqrt(sum((ua1-u4DRP1(end,:)).^2));

%Case 2
x2 = -200:dx:300;
t2 = [0 0 0 0:dt:200];

u2CDS2 = zeros(length(t2),length(x2)); %2nd Order Central Difference
u2CDS2(1,:) = x2>=-50 & x2<=50;
u2CDS2(2,:) = u2CDS2(1,:);
u2CDS2(3,:) = u2CDS2(1,:);
u2CDS2(4,:) = u2CDS2(1,:);
u4CDS2 = u2CDS2; %4th Order Central Difference
u6CDS2 = u2CDS2; %6th Order Central Difference
u4DRP2 = u2CDS2; %4th Order DRP
ua2 = double(x2>=150 & x2<=250); %analytic solution

for j = 4:length(t2)-1
    u2CDS2(j+1,:) = u2CDS2(j,:)-dt*b*NumericalDerivative(1,2,dx,u2CDS2(j:-1:j-3,:)); 
    u2CDS2(j+1,1) = 0;   %suppress left boundary fluctuation
    
    u4CDS2(j+1,:) = u4CDS2(j,:)-dt*b*NumericalDerivative(1,4,dx,u4CDS2(j:-1:j-3,:));
    u4CDS2(j+1,1:2) = 0;   %suppress left boundary fluctuation
    
    u6CDS2(j+1,:) = u6CDS2(j,:)-dt*b*NumericalDerivative(1,6,dx,u6CDS2(j:-1:j-3,:));
    u6CDS2(j+1,1:3) = 0;   %suppress left boundary fluctuation
    
    u4DRP2(j+1,:) = u4DRP2(j,:)-dt*b*FakeNumericalDerivative(dx,u4DRP2(j:-1:j-3,:));
    u4DRP2(j+1,1:3) = 0;   %suppress left boundary fluctuation    
end

%calc rms error
u2CDS2error = sqrt(sum((ua2-u2CDS2(end,:)).^2));
u4CDS2error = sqrt(sum((ua2-u4CDS2(end,:)).^2));
u6CDS2error = sqrt(sum((ua2-u6CDS2(end,:)).^2));
u4DRP2error = sqrt(sum((ua2-u4DRP2(end,:)).^2));

compute_time = toc;
save test2