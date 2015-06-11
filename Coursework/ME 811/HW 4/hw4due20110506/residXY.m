function [Ru,Rv,Riu,Riv] = residXY(U,AxO,AxE,AxW,AxN,AxS,Sx,V,AyO,AyE,AyW,AyN,AyS,Sy)

%this function calculates the residual for the X and Y momentum equations
%when given the link coefficients and the current velocity fields

tempU = zeros(size(U,1)+2,size(U,2)+2);
tempU(2:end-1,2:end-1) = U;
Riu = 0*U;
for i = 1:size(U,2)
	for j = 1:size(U,1)
		Riu(j,i) = Sx(j,i) - (AxO(j,i)*tempU(j+1,i+1) + ...
							  AxE(j,i)*tempU(j+1,i+2) + AxW(j,i)*tempU(j+1,i) + ...
							  AxN(j,i)*tempU(j+2,i+1) + AxS(j,i)*tempU(j,i+1));
	end
end
Ru = sqrt(sum(sum(Riu.^2)));

tempV = zeros(size(V,1)+2,size(V,2)+2);
tempV(2:end-1,2:end-1) = V;
Riv = 0*U;
for i = 1:size(V,2)
	for j = 1:size(V,1)
		Riv(j,i) = Sy(j,i) - (AyO(j,i)*tempV(j+1,i+1) + ...
							  AyE(j,i)*tempV(j+1,i+2) + AyW(j,i)*tempV(j+1,i) + ...
							  AyN(j,i)*tempV(j+2,i+1) + AyS(j,i)*tempV(j,i+1));
	end
end
Rv = sqrt(sum(sum(Riv.^2)));