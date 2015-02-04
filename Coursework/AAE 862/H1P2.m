function [] = H1P2()
%AAE 862 HW #1 Problem #1
%Completed by Michael Crawley
sigma = 0:.01:.9;
pratio = zeros(1,length(sigma));

for i = 1:length(sigma)
    pratio(i) = EQN(sigma(i));
end

hold on
plot(sigma,pratio);title('Pumping ratio versus Area ratio');xlabel('(Area Ratio)^-^1');ylabel('Pumping Ratio');
hold off

end

function [x] = EQN(sigma)
    a = (sigma/(1-sigma))^2+1;
    b = 4;
    c = -2*((1-sigma)/sigma);
    x = (-b + sqrt(b^2 -4*a*c))/(2*a);
end