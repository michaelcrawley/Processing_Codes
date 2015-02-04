function [Mach]=ShockJumpMach(M)
    gamma=1.4;
    Mach = sqrt((1+((gamma-1)/2)*M^2)/(gamma*M^2-(gamma-1)/2));
end