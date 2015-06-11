function [pp Resp] = solvePressEQ(Ap)
    Rtol = 1E-20;
    itrmax = 20;
    
    [pp ~] = ADIpm(Ap,Rtol,itrmax);
    [~, Resp] = calcRes(Ap,pp);
end