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