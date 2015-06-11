function [beta]=BetaSolver(M,theta)
    error = 1;
    i=1;
    alpha=45/2;
    beta = 45;
    while abs(error)>1E-100 && i<2000
        error = 2*cot(beta*pi/180)*((M^2*sin(beta*pi/180)^2-1)/(M^2*(1.4+cos(beta*pi/90))+2))-tan(theta*pi/180);
        beta = beta-sign(error)*alpha/i;
        i = i+1;
    end
end