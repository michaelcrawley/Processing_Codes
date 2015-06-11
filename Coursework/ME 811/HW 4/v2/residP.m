function [Rp,Rip] = residP(P,ApO,ApE,ApW,ApN,ApS,Sp)

%this function calculates the residual for the pressure correction equation
%when given the link coefficients and the current pressure field

tempP = zeros(size(P,1)+2,size(P,2)+2);
tempP(2:end-1,2:end-1) = P;
Rip = 0*P;
for i = 1:size(P,2)
    for j = 1:size(P,1)
        Rip(j,i) = Sp(j,i) - (ApO(j,i)*tempP(j+1,i+1) + ...
                              ApE(j,i)*tempP(j+1,i+2) + ApW(j,i)*tempP(j+1,i) + ...
                              ApN(j,i)*tempP(j+2,i+1) + ApS(j,i)*tempP(j,i+1));
    end
end
Rp = sqrt(sum(sum(Rip.^2)));