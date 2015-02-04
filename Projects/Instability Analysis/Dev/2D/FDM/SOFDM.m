function [phi,psi] = SOFDM(A,B,C,D,E,F,G,H,J,K,f,g,sigma,lambda,h)
    %A(phi)'' + B(phi)' + C(phi) + D(psi)' + E(psi) = f    
    %F(psi)'' + G(psi)' + H(psi) + J(phi)' + K(phi) = g    
    l = length(A);
    
    Z(1:2:2*l-1) = f;
    Z(2:2:2*l) = g;
    Z(1) = Z(1) - (A(1)/(h^2) - B(1)/(h*2))*sigma(1) + D(1)/(2*h)*lambda(1);
    Z(2) = Z(2) - (F(1)/(h^2) - G(1)/(h*2))*lambda(1) + J(1)/(2*h)*sigma(1);
    Z(end-1) = Z(end-1) -(A(end)/(h^2)+B(end)/(h*2))*sigma(2) - D(end)/(2*h)*lambda(2);
    Z(end) = Z(end) -(F(end)/(h^2)+G(end)/(h*2))*lambda(2) - J(end)/(2*h)*sigma(2);
    
    X1(1:2:2*l) = circshift(-J/(2*h),[0 -1]);
    X1(2:2:2*l) = 0;
    X2(1:2:2*l-1) = A/(h^2)-B/(2*h);
    X2(2:2:2*l) = F/(h^2)-G/(2*h); X2 = circshift(X2,[0 -2]);
    X3(1:2:2*l-1) = K;
    X3(2:2:2*l) = circshift(-D/(2*h),[0 -1]);
    X4(1:2:2*l-1) = C-2*A/(h^2);
    X4(2:2:2*l) = H-2*F/(h^2);
    X5(1:2:2*l-1) = E;
    X5(2:2:2*l) = J/(2*h); X5 = circshift(X5,[0 1]);
    X6(1:2:2*l-1)= A/(h^2)+B/(2*h);
    X6(2:2:2*l) = F/(h^2)+G/(2*h); X6 = circshift(X6,[0 2]);
    X7(1:2:2*l-1) = D/(2*h);
    X7(2:2:2*l) = 0; X7 = circshift(X7,[0 3]);
    X = spdiags([X1.' X2.' X3.' X4.' X5.' X6.' X7.'], -3:3,2*l,2*l);
    
    Y = X\(Z.');
    phi = Y(1:2:end-1);
    psi = Y(2:2:end);
    
end
