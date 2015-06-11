function [alpha vh] = ShearLayerLST(omega,N)
    y = linspace(-4,4,N);
    ly = length(y);
    dy = mean(diff(y));
    
    U = diag(0.5*(1+tanh(y)));
    U2 = diag(tanh(y).*(-(sech(y)).^2));
    
    d2vdy2 = mNumericalDerivative2(2,2,dy,ly);
    
    C0 = omega*d2vdy2;
    C1 = U.*d2vdy2-U2;
    C2 = diag(omega*ones(1,ly));
    C3 = -U;
    
    %solve Av = alpha*Mv
    A = [[-C1 -C2 -C3]; [eye(ly) zeros(ly) zeros(ly)]; [zeros(ly) eye(ly) zeros(ly)]];
    M = [[C0 zeros(ly) zeros(ly)]; [zeros(ly) eye(ly) zeros(ly)]; [zeros(ly) zeros(ly) eye(ly)]];

end