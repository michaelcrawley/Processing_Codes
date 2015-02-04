function [uh Resx]  = solveMomEQ(u,A,alpha)
    Rtol = 1E-20;
    itrmax = 2;
    [Ri Resx] = calcRes(A,u);    
    
    A.O = (1+alpha)*A.O;
    A.P = Ri;
    [up ~] = ADIpm(A,Rtol,itrmax);
    
    uh = u+up;
end