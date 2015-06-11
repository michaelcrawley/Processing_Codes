function [y,phi] = FDM(omg,alph,L,W,y1,y2)
    %Computations are performed using a second order finite difference
    %method
    h = .001; %Define step size
    y = y1:h:y2;
    l = length(y);
    U = 1+L*tanh(y)-W*exp(-log(2)*y.^2); %Define velocity profile
    U2 = -2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2); %Define d2/dy2 of velocity profile

    %Calculate Phi along y - Still need to calculate dPhi/dAlpha
    X = -(alph^2+U2./(U-omg/alph));
    Y = -(2*alph-(omg*U2)./((alph*U-omg).^2)); %calculate dX/dAlpha
    A = (1/h^2)*spdiags([ones(l,1) ((h^2)*X.'-2) ones(l,1)],-1:1,l,l);
    F = zeros(l,1);
    F(1) = -(1/h^2)*exp(alph*y1);F(end) = -(1/h^2)*exp(alph*-y2);
    phi = A\F;
end