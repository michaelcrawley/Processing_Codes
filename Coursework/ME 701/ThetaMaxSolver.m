function [thetamax]=ThetaMaxSolver(M)
    epsilon = 1E-3;
    beta=45;
    i=2;
    while (Theta(M,beta+epsilon)-Theta(M,beta))/epsilon > 1E-100  && i <2000
        deriv = (Theta(M,beta+epsilon)-Theta(M,beta))/epsilon;
        beta = beta+sign(deriv)*45/i;
    end
    thetamax = Theta(M,beta);
end

function [theta]=Theta(M,beta)
    theta = 180*atan(2*cot(beta*pi/180)*((M^2*sin(beta*pi/180)^2-1)/(M^2*(1.4+cos(beta*pi/90))+2)))/pi;
end