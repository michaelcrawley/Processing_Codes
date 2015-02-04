function [Rhalf] = calcRhalf(x,y,U)
    %Calc r = R (1/2) along the streamwise direction for axisymmetric jet
    %using PIV data
    %Inputs: x, y, U
    
%     [M N] = size(U);
%     if (length(x) ~= M && length(x) ~= N) || (length(y) ~= M && length(y) ~= N)
%         error('Input matrix size mismatch');
%     elseif (length(x) ~= N && length(x) == M)
%         U = U';
%     end
    Uc = max(U);
    Rhalf = zeros(1,length(x));
    for i = 1:length(x)
        Ut = U(:,i)/Uc(i);
        loc = find(Ut < 0.5,1,'last');  
        Rhalf(i) = interp1([Ut(loc+1) Ut(loc)],[y(loc+1) y(loc)], 0.5);
    end
end