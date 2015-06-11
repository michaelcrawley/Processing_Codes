function [Omega] = timestepping(bj,M,N,Delta)


    %Equation 1 - indirect
    options = optimset('TolFun',1e-15,'MaxFunEvals',1000,'MaxIter',1000);
    dO = 0.0001;
    Omega_bar = (5*dO:dO:0.3)';
    soln = zeros(length(Omega_bar)+1,2);
    fval = zeros(length(Omega_bar),1);
    extflag = fval;

    soln(1,:) = [Omega_bar(1),-Omega_bar(1)/1000];
    for n = 1:length(Omega_bar)
        eqn = @(x) abs(1i*(exp(-1i*(x(1)+1i*x(2)))-1)/(bj(1) + bj(2)*exp(1i*(x(1)+1i*x(2))) + bj(3)*exp(1i*2*(x(1)+1i*x(2))) + bj(4)*exp(1i*3*(x(1)+1i*x(2)))) - Omega_bar(n));
        [soln(n+1,:),fval(n),extflag(n)] = fminsearch(eqn,soln(n,:),options);
    end

    soln = soln(2:end,:);
    Omega1 = soln(:,2);

    %Equation 2 - direct
    C = 1.75*(M + sqrt(2))*20*N/(Delta*log(10)*(1-M));
    Omega2 = -Omega_bar/C;

    %Find intersection
    I = find(diff(sign(Omega2-Omega1)),1,'last')+1;
    Omega = Omega_bar(I);
end