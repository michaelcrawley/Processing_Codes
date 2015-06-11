function [Tf] = P3R2()
    ehs = 0.25; %emissivity of heating surface
    es = 0.7; %emissivity of slab
    ew = 0.1; %emissivity of walls
    h = 10; %W/m2/K
    Tw = 350; 
    Ths = 2000;
    Ta = 300;
    k = 30; %W/m/K
    sigma = 5.67E-8;
    dy = .001;

    e = [ew*ones(1,10) ehs ew*ones(1,9) es*ones(1,20) ew*ones(1,10)];

    Lx = [ 0:10:190 0:10:190 zeros(1,5) 200*ones(1,5) ]/1000 ;
    Rx = [10:10:200 10:10:200 zeros(1,5) 200*ones(1,5)]/1000;
    Ty = [zeros(1,20) 50*ones(1,20) 10:10:50 10:10:50 ]/1000;
    Ly = [zeros(1,20) 50*ones(1,20) 0:10:40 0:10:40 ]/1000;
    dx = mean(sqrt((Rx-Lx).^2+(Ty-Ly).^2));

    F = zeros(length(Lx));
    A = F; B = F;
    % elements in F: 
    %               1:20    top of oven (left to right)
    %               21:40   bottom of oven (left to right)
    %               41:45   left side of oven (bottom to top)
    %               46:50   right side of oven (bottom to top)
    for i = 1:length(Lx)
        for j = 1:length(Lx)
            if i ~= j
                d1 = sqrt((Rx(j)-Lx(i))^2+(Ly(i)-Ty(j))^2);
                d2 = sqrt((Lx(j)-Rx(i))^2+(Ly(j)-Ty(i))^2);
                s1 = sqrt((Rx(j)-Rx(i))^2+(Ty(j)-Ty(i))^2);
                s2 = sqrt((Lx(j)-Lx(i))^2+(Ly(j)-Ly(i))^2);
                F(i,j) = abs(d1+d2-s1-s2)/(2*dx);
            end
        end
    end
    
    Ti = mean([Tw Ths])*ones(20,1);
    Tf = fminsearch(@paramsearch,Ti);    
    function [error] = paramsearch(params)
        Tn = calcTn(F,e,sigma,Tw,Ths,Ta,params,dx,dy,h,k);
        error = sum(abs(Tn-params));
    end
end

function [Tn] = calcTn(F,e,sigma,Tw,Ths,Ta,Ti,dx,dy,h,k)
    [M N] = size(F);
    T = [Tw*ones(10,1); Ths; Tw*ones(9,1); Ti; Tw*ones(10,1)];
    
    A = zeros(size(F)); B = A;
    for i = 1:N
       for j = 1:M
            A(i,j) = KrD(i,j)/e(j)-(1/e(j)-1)*F(i,j);
            B(i,j) = KrD(i,j)-F(i,j);         
       end
    end
    Eb = sigma*T.^4;
    q = A\(B*Eb);
    D0 = [-h*(dx+dy)-k*dy/dx;(-2*k*dy/dx-h*dx)*ones(18,1);-h*(dx+dy)-k*dy/dx];
    D1 = k*dy/dx*ones(20,1);
    X = spdiags([D1 D0 D1],-1:1,20,20);
    Y = q(21:40)*dx-h*dx*Ta;
    Y(1) = Y(1)-h*dy*Ta;Y(end) = Y(end)-h*dy*Ta;
    Tn = X\Y;    
end
