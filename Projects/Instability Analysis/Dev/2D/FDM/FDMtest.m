function [] = FDMtest() 
    %Phi will be defined as sin(ßy)
    %Psi is then ycos(ßy)
    %and the ODEs are (phi)'' + ß^2(phi) = 0
    %and (psi)'' + ß^2(psi) + 2ß(phi) = 0 

    y1 = -pi; 
    y2 = pi;
    h = pi*.00001;
    y = y1:h:y2;
    l = length(y);
    
    f = y;
    g = zeros(1,l);
    
    A = ones(1,l);
    B = 2*ones(1,l);
    C = ones(1,l);
    D = zeros(1,l);
    E = D;
    F = D;
    G = D;
    H = D;
    J = D;
    K = D;
    
    sigma = [(exp(1-y1)+y1-2) (exp(1-y2)+y2-2)];
    lambda = [0 0];
    
    [phi,psi] = SOFDM(A,B,C,D,E,F,G,H,J,K,f,g,sigma,lambda,h);
    figure
    plot(y,phi,y,exp(1-y)+y-2);title('\Phi');legend('Numerical','True');
%     figure
%     plot(y,psi,y,0));title('\Psi');legend('Numerical','True');
end