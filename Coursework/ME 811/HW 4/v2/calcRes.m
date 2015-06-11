function [R2 Ri] = calcRes(Ao,Ae,Aw,An,As,S,phi)
%     Ae = [Ae(:,2:N) zeros(M,1)];
%     Aw = [zeros(M,1) Aw(:,1:N-1)];
%     An = [An(2:M,:); zeros(1,N)];
%     As = [zeros(1,N); As(1:M-1,:)];
%     
%     Ri = Ao.*phi+Ae.*phi+Aw.*phi+An.*phi+As.*phi-S;
    [M,N] = size(phi);
	Ri = As.*[zeros(1,N); phi(1:M-1,:)] ...		% phi(i,j-1)
		+ Aw.*[zeros(M,1), phi(:,1:N-1)] ...		% phi(i-1,j)
		+ Ao.*phi ...							% phi(i,j)
		+ Ae.*[phi(:,2:N), zeros(M,1)] ...		% phi(i+1,j)
		+ An.*[phi(2:M,:); zeros(1,N)] ...		% phi(i,j+1)
		- S;									% S(i,j)
     R2 = norm(Ri);
end