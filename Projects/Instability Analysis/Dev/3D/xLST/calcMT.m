function [theta] = calcMT(y,U)
    %Calculates momentum thickness based on incompressible flow
    theta = (U(2:end)'/Uc(1)).*(1-U(2:end)'/Uc(1))*abs(diff(y)');
end