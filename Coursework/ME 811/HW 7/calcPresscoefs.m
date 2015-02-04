function [Ap] = calcPresscoefs(uf,vf,rhof,dx,dy,xstepi,ystepi,Ao)
    [M N] = size(Ao);

    mimb = -(rhof.e.*uf.e-rhof.w.*uf.w)*dy-(rhof.n.*vf.n-rhof.s.*vf.s)*dx;
    
    Ap.O = zeros(M,N);
    Ap.E = Ap.O;
    Ap.W = Ap.O;
    Ap.N = Ap.O;
    Ap.S = Ap.O;
    Ap.P = -mimb;
    
    
    %interior nodes
    Ap.E(:,1:end-1) = 0.5*rhof.e(:,1:end-1)*dy*dy.*(1./Ao(:,1:end-1)+1./Ao(:,2:end));
    Ap.W(:,2:end) = 0.5*rhof.w(:,2:end)*dy*dy.*(1./Ao(:,2:end)+1./Ao(:,1:end-1));
    Ap.N(1:end-1,:) = 0.5*rhof.n(1:end-1,:)*dx*dx.*(1./Ao(1:end-1,:)+1./Ao(2:end,:));
    Ap.S(2:end,:) = 0.5*rhof.s(2:end,:)*dx*dx.*(1./Ao(2:end,:)+1./Ao(1:end-1,:));
    
    %step wall nodes
    Ap.W(1:ystepi,xstepi+1) = 0;
    Ap.S(ystepi+1,1:xstepi) = 0;
    
    Ap.O = -(Ap.E+Ap.W+Ap.N+Ap.S);
    
    %outlet nodes
    Ap.O(:,end) = Ap.O(:,end)-rhof.e(:,end)*dy*dy.*(1./Ao(:,end));
    
    %In step nodes
    Ap.O(1:ystepi,1:xstepi) = 1;
    Ap.E(1:ystepi,1:xstepi) = 0;
    Ap.W(1:ystepi,1:xstepi) = 0;
    Ap.N(1:ystepi,1:xstepi) = 0;
    Ap.S(1:ystepi,1:xstepi) = 0;
    Ap.P(1:ystepi,1:xstepi) = 0;
end