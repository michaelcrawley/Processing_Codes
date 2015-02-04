function [M2]=ObliqueShockJumpMach(M,theta,beta)
    M1N = M*sin(beta*pi/180);
    M2N = ShockJumpMach(M1N);
    M2 = M2N/sin((beta-theta)*pi/180);
end