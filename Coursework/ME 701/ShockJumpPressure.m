function [pratio]=ShockJumpPressure(M)
    gamma=1.4;
    pratio=1+((2*gamma)/(gamma+1))*(M^2-1);
end