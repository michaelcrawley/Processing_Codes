function [pratio] = testfunction(M)
gamma = 1.4;
    if M >= 1
        pratio = ((((gamma+1)*M^2)/(2+(gamma-1)*M^2))^(gamma/(gamma-1)))/((1+2*gamma/(gamma+1)*(M^2-1))^(1/(gamma-1))) * (1+0.5*(gamma-1)*M^2)^(gamma/(gamma-1));         
    else
        pratio = (1+0.5*(gamma-1)*M^2)^(gamma/(gamma-1));
    end

end