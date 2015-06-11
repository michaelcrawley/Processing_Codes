function [val] = viscous(T1,T2,u1,u2,x,dx)
    muo = 1.716E-5; %kg/ms
    To = 273; %K
    T = T1+(T2-T1)*(dx./x);
    dudx = (u2-u1)/x;    
    val = (4/3)*muo*dudx^2*(1./T).*(T/To).^0.77;
end