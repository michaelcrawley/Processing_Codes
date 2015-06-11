function []=P1(M1,theta2,theta3, P1)
    epsilon = 1E-3;
    beta2 = BetaSolver(M1,theta2);
    beta3 = BetaSolver(M1,theta3);
    [M2 M12N]=ObliqueShockJumpMach(M1,theta2, beta2);
    [M3 M13N]=ObliqueShockJumpMach(M1,theta3, beta3);
    Pr2 = ShockJumpPressure(M12N);
    Pr3 = ShockJumpPressure(M13N);    
    Phie = 0;
    error=1;
    counter = 0;
    while abs(error)>1E-6 || counter < 1000
        [error Pr4 Pr4p] = Phi(M2,M3,Pr2,Pr3,theta2,theta3,Phie);
        errorrate = (Phi(M2,M3,Pr2,Pr3,theta2,theta3,Phie+epsilon)-error)/epsilon;
        Phie = Phie-error/errorrate;
        counter = counter +1;
    end
    P4 = Pr3*Pr4*P1;
    P4p = Pr2*Pr4p*P1;
    fprintf('The Angle of the Slipline is: %6.5f degrees \n',Phie);
    fprintf('The Pressure in Region 4 is: %6.5f atm \n',P4);
    fprintf('The Pressure in Region 4 prime is: %6.5f atm \n',P4p);
end

function [M2, M1N]=ObliqueShockJumpMach(M,theta,beta)
    M1N = M*sin(beta*pi/180);
    M2N = ShockJumpMach(M1N);
    M2 = M2N/sin((beta-theta)*pi/180);
end

function [pratio]=ShockJumpPressure(M)
    gamma=1.4;
    pratio=1+((2*gamma)/(gamma+1))*(M^2-1);
end

function [Mach]=ShockJumpMach(M)
    gamma=1.4;
    Mach = sqrt((1+((gamma-1)/2)*M^2)/(gamma*M^2-(gamma-1)/2));
end

function [error, Pr4, Pr4p]=Phi(M2,M3,Pr2,Pr3,theta2,theta3,Phi)
    theta4 = theta3+Phi;
    theta4p= theta2-Phi;
    beta4 = BetaSolver(M3,theta4);
    beta4p = BetaSolver(M2,theta4p);
    [M4 M3N]=ObliqueShockJumpMach(M3,theta4,beta4);
    [M4p M2N]=ObliqueShockJumpMach(M2,theta4p,beta4p);
    Pr4 = ShockJumpPressure(M3N);
    Pr4p = ShockJumpPressure(M2N);
    error = (Pr4*Pr3)-(Pr4p*Pr2);
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