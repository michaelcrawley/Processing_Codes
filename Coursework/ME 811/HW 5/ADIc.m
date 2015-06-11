function [phi R2] = ADIc( A,B,C,D,E,F , phi , tol , maxIt )
% This function employs the ADI method to solve a system of equations on a
% non-uniform 2D mesh grid. Dirichlet boundary conditions are assumed and
% must be enforced in the pentadiagonal matrix.
%	Inputs are of the form:
%		A	[1 K]or[M N]	-N sub-diagonal
%		B	[1 K]or[M N]	-1 sub-diagonal
%		C	[1 K]or[M N]	Main-diagonal
%		D	[1 K]or[M N]	+1 super-diagonal
%		E	[1 K]or[M N]	+N super-diagonal
%		F	[1 K]or[M N]	Source term
%		phi		[M N]		Initial guess
%		tol		[1 1]		Absolute (+), or relative (-) tolerance

%% Determine the 2D grid dimensions
[M N] = size(phi);
K = M*N;

%% Check for proper inputs
if ~( ...
		numel(A)==K && numel(B)==K && ...
		numel(C)==K && numel(D)==K && ...
		numel(E)==K && numel(F)==K ...
	)
	fprintf(2,'Inputs do not match... terminating.\n');
	R2 = -1;
	return;
end

%% Rearrange pentadiagonal matrix into grid form
if size(A)~=[M N], A = reshape(A,M,N); end
if size(B)~=[M N], B = reshape(B,M,N); end
if size(C)~=[M N], C = reshape(C,M,N); end
if size(D)~=[M N], D = reshape(D,M,N); end
if size(E)~=[M N], E = reshape(E,M,N); end
if size(F)~=[M N], F = reshape(F,M,N); end

%% Clear non-existant cells
A(1,:) = 0;		B(:,1) = 0;
D(:,N) = 0;		E(M,:) = 0;

%% Establish conditions for convergence
if ~exist('tol','var')
	% Reduce residual by two orders of magnitude as a default
	tol = -2;
end
if tol<=0
	% Determine cutoff for relative convergence
	Ri = A.*[zeros(1,N); phi(1:M-1,:)] ...		% phi(i,j-1)
		+ B.*[zeros(M,1), phi(:,1:N-1)] ...		% phi(i-1,j)
		+ C.*phi ...							% phi(i,j)
		+ D.*[phi(:,2:N), zeros(M,1)] ...		% phi(i+1,j)
		+ E.*[phi(2:M,:); zeros(1,N)] ...		% phi(i,j+1)
		- F;									% S(i,j)
	R2 = sqrt(sum(sum( Ri.^2 )));
	tol = R2*10^tol;
end

%% Begin outer loop
if ~exist('maxIt','var'), maxIt=1000; end
for i=1:maxIt
	if rem(i,2)==1
		% Perform row-wise sweep
		S = F(1,:) - E(1,:).*phi(2,:);
		phi(1,:) = TDMA(B(1,:),C(1,:),D(1,:),S);
		for m=2:M-1
			S = F(m,:) - A(m,:).*phi(m-1,:) - E(m,:).*phi(m+1,:);
			phi(m,:) = TDMA(B(m,:),C(m,:),D(m,:),S);
		end
		S = F(M,:) - phi(M-1,:).*A(M,:);
		phi(M,:) = TDMA(B(M,:),C(M,:),D(M,:),S);
	else
		% Perform column-wise sweep
		S = F(:,1) - D(:,1).*phi(:,2);
		phi(:,1) = TDMA(A(:,1),C(:,1),E(:,1),S);
		for n=2:N-1
			S = F(:,n) - B(:,n).*phi(:,n-1) - D(:,n).*phi(:,n+1);
			phi(:,n) = TDMA(A(:,n),C(:,n),E(:,n),S);
		end
		S = F(:,N) - B(:,N).*phi(:,N-1);
		phi(:,N) = TDMA(A(:,N),C(:,N),E(:,N),S);
	end
	
	% Calculate the residual
	Ri = A.*[zeros(1,N); phi(1:M-1,:)] ...		% phi(i,j-1)
		+ B.*[zeros(M,1), phi(:,1:N-1)] ...		% phi(i-1,j)
		+ C.*phi ...							% phi(i,j)
		+ D.*[phi(:,2:N), zeros(M,1)] ...		% phi(i+1,j)
		+ E.*[phi(2:M,:); zeros(1,N)] ...		% phi(i,j+1)
		- F;									% S(i,j)
	R2(i) = sqrt(sum(sum( Ri.^2 )));
	
	% Check for convergence
	if R2(i)<=tol, return; end
end

% % % % If the loop completes, then the solution did not converge
% % % fprintf(2,'Solution did not converge.\n');
