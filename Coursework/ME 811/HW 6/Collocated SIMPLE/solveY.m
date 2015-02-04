function [v R2] = solveY(A,v,al,iter)
% Solve the x-momentum equation in correction form using the ADI method.

[M N] = size(A.o);

% Calculate the residual
Ri = A.s .* [zeros(1,N); v(1:M-1,:)] ...		% v(i,j-1)
	+ A.w .* [zeros(M,1), v(:,1:N-1)] ...		% v(i-1,j)
	+ A.o.* v ...								% v(i,j)
	+ A.e .* [v(:,2:N), zeros(M,1)] ...			% v(i+1,j)
	+ A.n .* [v(2:M,:); zeros(1,N)] ...			% v(i,j+1)
	- A.p.y;									% S(i,j)
R2 = sqrt(sum(sum( Ri.^2 )));

% Initial guess
null = zeros(M,N);

% Solve
dv = ADI( A.s , A.w , (1+al)*A.o , A.e , A.n , -Ri , null , -inf , iter );

% Update
v = v + dv;