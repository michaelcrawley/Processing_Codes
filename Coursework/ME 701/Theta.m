function [theta]=Theta(M,beta)
    theta = 180*atan(2*cot(beta*pi/180)*((M^2*sin(beta*pi/180)^2-1)/(M^2*(1.4+cos(beta*pi/90))+2)))/pi;
end