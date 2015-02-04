function []=P2(M,dSN)
    epsilon = .1;
    SN = dSN+1;
    i=1;
    while SN >= dSN 
        SN = CalcShockNumber(M,i*epsilon);
        i = i+1;
    end
    maxdeflection = (i-2)*epsilon;
    fprintf('Maximum deflection angle for which %i oblique shocks occur is: %2.1f degrees\n', dSN, maxdeflection);
end

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

function [SN]=CalcShockNumber(M,theta)
    SN = 0;
    thetamax = ThetaMaxSolver(M);
    while theta<thetamax
        beta = BetaSolver(M,theta);
        M = ObliqueShockJumpMach(M,theta,beta);
        thetamax = ThetaMaxSolver(M);
        SN = SN+1;
    end
end

function [M2]=ObliqueShockJumpMach(M,theta,beta)
    M1N = M*sin(beta*pi/180);
    M2N = ShockJumpMach(M1N);
    M2 = M2N/sin((beta-theta)*pi/180);
end

function [Mach]=ShockJumpMach(M)
    gamma=1.4;
    Mach = sqrt((1+((gamma-1)/2)*M^2)/(gamma*M^2-(gamma-1)/2));
end

function [beta]=BetaSolver(M,theta)
    error = 1;
    i=1;
    alpha=45;
    beta = 1;
    while abs(error)>1E-100 && i<1000
        error = 2*cot(beta*pi/180)*((M^2*sin(beta*pi/180)^2-1)/(M^2*(1.4+cos(beta*pi/90))+2))-tan(theta*pi/180);
        beta = beta-sign(error)*alpha/i;
        i = i+1;
    end
end