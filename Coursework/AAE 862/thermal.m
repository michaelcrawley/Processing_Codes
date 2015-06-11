function [val] = thermal(T1,T2,x,dx)
    muo = 1.716E-5; %kg/ms
    To = 273; %K
    cp = 1004.5;
    Pr = 0.71;
    T = T1+(T2-T1)*(dx./x);
    dTdx = (T2-T1)/x;
    val = (cp/Pr)*muo*dTdx^2*(1./(T.^2)).*(T/To).^0.77;
end