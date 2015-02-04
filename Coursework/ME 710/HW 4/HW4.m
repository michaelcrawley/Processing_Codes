function [] = HW4(xs)
%Completed by Michael Crawley for ME 710 HW#4
    etamax = 7;
    Pr = 6;
    Re = 400;    
    eta = 0:0.01:etamax;
    f = @(x) (0.75*(x.^2/4.64) - 0.125*(x.^4/(4.64^3))).*(x<=4.64)+(2.9+x-4.64).*(x > 4.64);
    p = @(x,Pr) exp(-0.5*Pr*cumtrapz(x,f(x)));
    
   
    Po = -1/(trapz(eta,p(eta,Pr)));
    theta =1+Po*cumtrapz(eta,p(eta,Pr)); %(T-Tinf)/(Tw-Tinf)
    ys = zeros(length(xs),length(eta));
    figure; hold on;
    for i = 1:length(xs)
        ys(i,:) = eta/sqrt(Re/xs(i));
        plot(theta,ys(i,:),'Color',rand(1,3));
    end
    ylabel('y*');title('Nondimensionalized Temperature Distribution');xlabel('\theta');legend(num2str(xs'));

    xs = 0.001:0.001:1;
    Nux = zeros(1,length(xs));
    for i = 1:length(xs)
        yst = eta/sqrt(Re/xs(i));
        dtheta = NumericalDerivative(1,1,mean(diff(yst)),theta);
        Nux(i) = -xs(i)*dtheta(1)/theta(1);
    end
    figure;
    plot(xs,Nux);ylabel('Nu_x');xlabel('x/L');title('Local Nusselt number along the flat plat');
end
