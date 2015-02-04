function [T,Z] = H1P1()    
    zo = [ 0.1 0.1 0.1];
    [T,Z] = ode113(@Lorenz,[0 100], zo);
    figure;
    subplot(3,1,1), plot(T,Z(:,1)); xlabel('t'); ylabel('x');
    subplot(3,1,2), plot(T,Z(:,2)); xlabel('t'); ylabel('y');
    subplot(3,1,3), plot(T,Z(:,3)); xlabel('t'); ylabel('z');
    figure;
    plot3(Z(:,1),Z(:,2),Z(:,3));

end

function [dz] = Lorenz(t,z)
 	rho = 28;
    sigma = 10; beta = 8/3;
    dz = zeros(3,1);
    dz(1) = sigma*(z(2)-z(1));
    dz(2) = z(1)*rho-z(2)-z(1)*z(3);
    dz(3) = z(1)*z(2)-beta*z(3);

end