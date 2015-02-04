function [y,phi,A,F,X] = temp_FDM(omg,alph,L,W,y1,y2)
    h = .001; %Define step size
    y = y1:h:y2;
    l = length(y);
    U = 1+L*tanh(y)-W*exp(-log(2)*y.^2); %Define velocity profile
    U2 = -2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2); %Define second derivative of velocity profile
    X = -(alph^2+U2./(U-omg/alph));
    A = (1/h^2)*spdiags([ones(l,1) (h^2)*X'-2 ones(l,1)],-1:1,l,l);
    F = zeros(l,1);
    F(1) = -(1/h^2)*exp(alph*y1);F(end) = -(1/h^2)*exp(alph*-y2);
    phi = A\F;
end