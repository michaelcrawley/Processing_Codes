function [] = H1P1()
%AAE 862 HW #1 Problem #1
%Completed by Michael Crawley

T = 234.7; %K
gamma = 1.4;
c = sqrt(gamma*287*T);
M = 0:.01:2;
C_pratio = zeros(1,length(M));
C_pratio2 = zeros(1,length(M));

for i = 1:length(M)
    C_pratio(i) = Compressible(M(i));
end
X1 = find(C_pratio > 3, 1 );
for i = 1:length(M)
    C_pratio2(i) = Bernoulli(M(i));
end
X2 = find(C_pratio2 > 3, 1 );

plot(C_pratio(1:X1),M(1:X1), C_pratio2(1:X2), M(1:X2));xlim([1 3]);title('Comparison of Bernoulli Eqn versus Actual Measured Velocity');
xlabel('Pt/P');ylabel('Mach Number');legend('Measured', 'Bernoulli');

end

function [pratio] = Compressible(M)
    gamma = 1.4;
    if M >= 1
        pratio = ((((gamma+1)*M^2)/(2+(gamma-1)*M^2))^(gamma/(gamma-1)))/((1+2*gamma/(gamma+1)*(M^2-1))^(1/(gamma-1))) * (1+0.5*(gamma-1)*M^2)^(gamma/(gamma-1));         
    else
        pratio = (1+0.5*(gamma-1)*M^2)^(gamma/(gamma-1));
    end
end

function [pratio] = Bernoulli(M)
    gamma = 1.4;
    R = 287;
    T = 234.7;
    c = sqrt(gamma*287*T);
    pratio = 1+0.5*(1/(R*T))*(M*c)^2;
end