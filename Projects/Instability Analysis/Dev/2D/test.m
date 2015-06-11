function dz=test(y,z,alpha,c)
    % y - independent variable
    % z - dependent variable
    % alpha - complex wave number
    % c = omega/alpha
    % dz = dz/dy at y
    L=0.45;
    W=0.8;
    u=1+L*tanh(y)-W*exp(-log(2)*y.^2);
    u2=-2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2);
    dz=zeros(4,1);
    %[u,u2]=u_sub(y);
    cmc=u-c;
    dz(1)=z(2);                        
    dz(2)=(alpha*alpha+u2/cmc)*z(1) ;  
    dz(3)=z(4);						   
    dz(4)=(alpha*alpha+u2/cmc)*z(3)+z(1)*(2*alpha-c*u2/(cmc*cmc*alpha));
end