function [ds] = NShock(M)
    cp = 1004.5;
    R = 287;
    Tr = (1+2*(1.4/2.4)*(M^2-1))*((2+0.4*M^2)/(2.4*M^2));
    Pr = 1+2*(1.4/2.4)*(M^2-1);
    ds = cp*log(Tr)-R*log(Pr);
end