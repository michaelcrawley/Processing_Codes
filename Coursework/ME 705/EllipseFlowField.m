function [u,v]=EllipseFlowField(U,a,c,x,y)
    z=complex(x,y);
    W=U*(1+((a/c)^2-1)*(0.5-z/(4*sqrt((z/2)^2-c^2))));
    u=real(W);
    v=-imag(W);
end