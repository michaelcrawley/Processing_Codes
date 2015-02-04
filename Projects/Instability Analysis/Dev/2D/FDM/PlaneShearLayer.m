function [phi,psi,y,X,Y,Z, C, K,er,der,dalph] = PlaneShearLayer(omg,alph,h)
    %Computations are performed using a second order finite difference
    %method
    L = 0.45;
    W = 0.8;

    y1 = -8; 
    y2 = 8;
    y = y1:h:y2;
    l = length(y);
    U = 1+L*tanh(y)-W*exp(-log(2)*y.^2); %Define velocity profile
    U2 = (-2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2)); %Define d2/dy2 of velocity profile
    
    A = ones(1,l);
    B = zeros(1,l);
    C = -(alph^2+U2./(U-omg/alph));
    D = zeros(1,l);
    E = zeros(1,l);
    F = ones(1,l);
    G = zeros(1,l);
    H = C;
    J = zeros(1,l);
    K = -(2*alph-(omg*U2)./((U*alph-omg).^2));
    
    f = zeros(1,l);
    g = zeros(1,l);
    
    sigma = [exp(y1*alph) exp(-y2*alph)];%Boundary Conditions for Phi
    lambda = [y1*exp(y1*alph) -y2*exp(-y2*alph)];%Boundary Conditions for Psi
    
    [phi,psi,X,Y,Z] = SOFDM(A,B,C,D,E,F,G,H,J,K,f,g,sigma,lambda,h);
    phi = phi.';
    psi = psi.';    

    er1 = -alph*phi(1)+1/(2*h)*(phi(2)-sigma(1));
    er2 = alph*phi(end)+1/(2*h)*(sigma(2)-phi(end-1));
    er = [ er1; er2];
    
    der1 = -phi(1)-alph*psi(1)+1/(2*h)*(psi(2)-lambda(1));
    der2 = phi(end)+alph*psi(end)+1/(2*h)*(lambda(2)-psi(end-1));
    der = [ der1; der2];
    
    dalph1=-er1/der1;
    dalph2=-er2/der2;
    dalph = [ dalph1; dalph2];
    
    figure
    plot(y,real(phi),y,imag(phi));title('\Phi');
    figure
    plot(y,real(psi),y,imag(psi));title('\Psi');
end