function [u R2] = solveX(A,u,al,iter)
% Solve the x-momentum equation in correction form using the ADI method.

[M N] = size(A.o);

% Calculate the residual
Ri = A.s .* [zeros(1,N); u(1:M-1,:)] ...		% u(i,j-1)
	+ A.w .* [zeros(M,1), u(:,1:N-1)] ...		% u(i-1,j)
	+ A.o.* u ...								% u(i,j)
	+ A.e .* [u(:,2:N), zeros(M,1)] ...			% u(i+1,j)
	+ A.n .* [u(2:M,:); zeros(1,N)] ...			% u(i,j+1)
	- A.p.x;									% S(i,j)
R2 = sqrt(sum(sum( Ri.^2 )));

% Initial guess
null = zeros(M,N);

% Solve
du = ADI( A.s , A.w , (1+al)*A.o , A.e , A.n , -Ri , null , -inf , iter );

% Update
u = u + du;