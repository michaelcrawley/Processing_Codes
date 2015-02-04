%=================================
function dz=deriv(y,z,alpha,c,m,nmode)
% y - independent variable
% z - dependent variable
% alpha - complex wave number
% c = omega/alpha
% dz = dz/dy at y
dz=zeros(4,1);
[u,u1]=u_sub(y);
[tep,tep1]=t_sub(y,m);
cmc=u-c;
dz(1)=z(2);                        
dz(2)=(2*u1/cmc-tep1/tep-1./y)*z(2)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(1);  
dz(3)=z(4);						   
dz(4)=(2*u1/cmc-tep1/tep-1./y)*z(4)+2*alpha*(1-m*m*cmc*cmc/tep)*z(1)+nmode^2/y^2*z(3)+alpha*alpha*(1-m*m*cmc*cmc/tep)*z(3)...
    -z(1)*2*m^2*cmc*c*alpha/tep-2*u1*c*z(2)/(cmc*cmc*alpha);




%=================================

