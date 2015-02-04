function [phi,R] = ADIn(A,B,C,D,E,S,guess,Rtol)

%ADI uses the "Alternating Direction Implicit" method to solve the
%specified set of linear algebraic equations when given the 5 diagonals of
%the pentadiagonal system, with "C" being the diagonal, "B", and "D" being
%the off diagonals, and "A" and "E" being the separate bands. The
%right-hand side vector is also necessary, along with the residual
%tolerance: a positive number is used as an absolute tolerance for the
%residual, while a negative number represents the number of orders of
%magnitude for the residual to be reduced. If the argument is a string, it
%represents the number of sweeps to be performed. The code returns the
%solution, and a vector of the residual values so that the system
%convergence may be studied. An example use is shown below:
%
%[soln,resid] = ADI(A,B,C,D,E,S,guess,Rtol)
%
%where guess is an initial guess
%note that this code assumes a counting pattern along columns first,
%starting with the top row of the matrix

%global phi R n

%parsing the inputs
N = size(guess,1);
M = size(guess,2);

K = M*N;
if size(C,1) <= size(C,2)
	A = A';
	B = B';
	C = C';
	D = D';
	E = E';
	S = S';
end

if length(S) ~= K || length(C) ~= K
	fprintf('Error: Vector length does not match given grid parameters.\n\n');
	return
end

if length(B) ~= K
	B = [0;B];
end
if length(D) ~= K
	D = [D;0];
end

if length(A) ~= K
	A = [zeros(M,1);A];
	A = A(1:K);
end
if length(E) ~= K
	E = [zeros(M,1);E];
	E = E(1:K);
end

A = reshape(A,M,N)';
B = reshape(B,M,N)';
C = reshape(C,M,N)';
D = reshape(D,M,N)';
E = reshape(E,M,N)';
S = reshape(S,M,N)';

%initializing the solution matrix and padding with zeros to ease coding,
%these will be cropped later
phi = guess;
phi = [zeros(1,M);phi;zeros(1,M)];
phi = [zeros(N+2,1),phi,zeros(N+2,1)];

%calculating the initial residual
Ri = A.*phi(1:N,2:M+1) + B.*phi(2:N+1,1:M) + C.*phi(2:N+1,2:M+1) + D.*phi(2:N+1,3:M+2) + E.*phi(3:N+2,2:M+1) - S;
R = sqrt(sum(sum(Ri.^2)));

%preprocessing the residual tolerance argument and determining an absolute
%tolerance
if ischar(Rtol)
    sweeps_lim = str2num(Rtol);
    Rtol = 0;
else
    if Rtol > 0
        Rtol = Rtol;
    else
        Rtol = R*10^Rtol;
    end
    sweeps_lim = inf;
end

sweeps = 0;
n = 2;
while R(n-1) >= Rtol && sweeps < sweeps_lim
	%sweeping to update phi
	for i = 1:N
		diag = C(i,1:M);
		sub = B(i,1:M);
		sup = D(i,1:M);
		RHS = S(i,1:M) - (A(i,1:M).*phi(i,2:M+1) + E(i,1:M).*phi(i+2,2:M+1));

		phi(i+1,2:M+1) = TDMS(diag,RHS,sub,sup)';

	end

	for j = 1:M
		diag = C(1:N,j);
		sub = A(1:N,j);
		sup = E(1:N,j);
		RHS = S(1:N,j) - (B(1:N,j).*phi(2:N+1,j) + D(1:N,j).*phi(2:N+1,j+2));

		phi(2:N+1,j+1) = TDMS(diag,RHS,sub,sup);

	end

	Ri = A.*phi(1:N,2:M+1) + B.*phi(2:N+1,1:M) + C.*phi(2:N+1,2:M+1) + D.*phi(2:N+1,3:M+2) + E.*phi(3:N+2,2:M+1) - S;
	R(n) = sqrt(sum(sum(Ri.^2)));
	
	n = n + 1;
    sweeps = sweeps + 1;
	
end

%postprocessing R
if size(R,1) <= size(R,2)
	R = R';
end

%cropping the solution back to its original size
phi = phi(2:N+1,2:M+1);