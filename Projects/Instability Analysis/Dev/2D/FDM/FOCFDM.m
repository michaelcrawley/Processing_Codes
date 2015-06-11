function [y,phi] = FOCFDM(alpha,omg,h)
    L = 0.45;
    W = 0.8;
    
    y1 = -8; 
    y2 = 8;
    y = y1:h:y2;
    l = length(y);
    
    U = 1+L*tanh(y)-W*exp(-log(2)*y.^2);
    U1 = L*(sech(y).^2)+W*log(2).*y.*2.^(1-(y.^2));
    U2 = -2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2);
    U3 = 8*W*log(2)^3.*y.^3.*2.^(-(y.^2))-12*W*log(2)^2.*y.*2.^(-(y.^2))+4*L.*sech(y).^2.*tanh(y).^2-2*L.*sech(y).^4;
    U4 = -16*W*log(2)^4.*y.^4.*2.^(-(y.^2))+48*W*log(2)^3.*y.^2.*2.^(-(y.^2))-12*W*log(2)^2.*2.^(-(y.^2))-8*L.*sech(y).^2.*tanh(y).^3+16*L*sech(y).^4.*tanh(y);

    X = -(alpha^2+U2./(U-omg/alpha));
    Y = -(U3./(U-omg/alpha)-(U1.*U2)./((U-omg/alpha).^2));
    Z = -(U4./(U-omg/alpha)-(2*U1.*U3)./((U-omg/alpha).^2)+(2*(U1.^2).*U2)./((U-omg/alpha).^3)-(U2./(U-omg/alpha)).^2);
    
    A = 1/h^2+1/12-1/(12*h)*Y;
    B = X +(h^2)*Z/12-2/(h^2)-1/6;
    C = 1/h^2+1/12+1/(12*h)*Y;

    f = zeros(l,1);
    sigma = [exp(y1*alpha) exp(-y2*alpha)];%Boundary Conditions for Phi
    f(1) = -A(1)*sigma(1);
    f(end) = -C(end)*sigma(end);
    
    U = spdiags([ circshift(A,[0 -1])' B' circshift(C,[0 1])'],-1:1,l,l);
    phi = U\f;
end