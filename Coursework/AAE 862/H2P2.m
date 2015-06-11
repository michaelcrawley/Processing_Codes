function [] = H2P2(M1)
    To = 273; %K    
    Te = @(To,M) To/(1+0.2*M^2);    
    
    %Part C 
    M2 = M_NShock(M1); 
    T1 = Te(To,M1); 
    T2 = Te(To,M2); 
    u1 = M1*sqrt(1.4*287*T1); 
    u2 = M2*sqrt(1.4*287*T2); 
    ds = S_NShock(M1); %Entropy generation from 1D shock relations 
    
    error = 1;
    dx = 1E-7;
    x = 0;
    n = 0;
    while(abs(error)>1E-3) && n<1E4
        x = x+dx;
        dsi = quad(@(dx) thermal(T1,T2,x,dx) ,0,x) + quad(@(dx) viscous(T1,T2,u1,u2,x,dx), 0, x); 
        error = (ds-dsi)/ds;
        n = n+1;
    end
    x
    x = 0:dx:x;
    s_t = zeros(1,length(x));
    s_v = zeros(1,length(x));
    for i = 1:length(x)
        s_t(i) = quad(@(dx) thermal(T1,T2,x(end),dx) ,0,x(i))/ds;
        s_v(i) = quad(@(dx) viscous(T1,T2,u1,u2,x(end),dx), 0, x(i))/ds;
    end
    s = s_t+s_v;
    figure; plot(x/max(x),s,x/max(x),s_t,x/max(x),s_v); xlim([0 1]); ylim([0 1]); legend('Total','Heat Transfer','Viscous Dissipation'); ylabel('(s-s1)/(s2-s1)'); xlabel('x/t');  
end


function [M2] = M_NShock(M)
    M2 = sqrt((1+((1.4-1)/2)*M^2)/(1.4*M^2-0.5*0.4));
end

function [ds] = S_NShock(M)
    cp = 1004.5;
    R = 287;
    Tr = (1+2*(1.4/2.4)*(M^2-1))*((2+0.4*M^2)/(2.4*M^2));
    Pr = 1+2*(1.4/2.4)*(M^2-1);
    ds = cp*log(Tr)-R*log(Pr); %in J/kg*K
end

function [val] = thermal(T1,T2,x,dx)
    muo = 1.716E-5; %kg/ms
    To = 273; %K
    cp = 1004.5;
    Pr = 0.71;
    T = T1+(T2-T1)*(dx./x);
    dTdx = (T2-T1)/x;
    val = (cp/Pr)*muo*dTdx^2*(1./(T.^2)).*(T/To).^0.77;
end

function [val] = viscous(T1,T2,u1,u2,x,dx)
    muo = 1.716E-5; %kg/ms
    To = 273; %K
    T = T1+(T2-T1)*(dx./x);
    dudx = (u2-u1)/x;    
    val = (4/3)*muo*dudx^2*(1./T).*(T/To).^0.77;
end