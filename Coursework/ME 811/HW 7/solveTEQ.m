function [Th ResT]  = solveTEQ(T,A,alpha)
    Rtol = 1E-20;
    itrmax = 20;
    [Ri ResT] = calcRes(A,T);    
    
    A.O = (1+alpha)*A.O;
    A.P = Ri;
    [Tp ~] = ADIpm(A,Rtol,itrmax);
    
    Th = T+Tp;
end