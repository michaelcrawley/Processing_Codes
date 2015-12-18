function [phi,R] = ADI2d(X,S,method,Rtol,itrmax)
    %Performs the Alternating Direction Implicit solver for a 2 dimensional
    %system.
    %Note:  conversion from 2-D grid to 1-D array must be done as i,j -> k = M(i-1)+j,
    %where there are M nodes in the y direction and N nodes in the x
    %direction
    %Code Version: 1.0 @ 2011-04-10
    %Inputs:
    %       X:      matrix for the x-derivative terms (M*N,M*N)
    %       S:      source term, either matrix or scalar (M,N)
    %       method: '-TDMA' for tridiag systems, '-GaussSeidel' otherwise
    %       Rtol:   Residual tolerance (optional; default 1E-5)
    %       itrmax: maximum number of iterations (optional; default 1E4)
    %Outpus:
    %       phi:    Solution to the problem
    %       R:      L2norm of the Residuals at each iteration
    
    if ~exist('Rtol','var')
        Rtol = 1E-5;
    end
    if ~exist('itrmax','var')
        itrmax = 1E4;
    end    
    
    [M N] = size(S);
    phi = zeros(M*N,1);
    S = reshape(S,M*N,1);
 
    R = zeros(1,itrmax);
    R(1) = norm(X*phi-S);
    counter = 0;
    while (R(counter+1) >= Rtol) && (counter < itrmax)        
        for j = 1:M %Row sweep
            k = j:M:j+M*(N-1); %determine nodal array points
            Yt = X; Yt(k,k) = 0;
            t = Yt*phi;
            b = S(k)-t(k);
            A = X(k,k);
            if strcmpi(method,'-TDMA') %use tridiagonal matrix solver
                phi(k) = TDMsolver(A,b);
            elseif strcmpi(method,'-GaussSeidel') %use Gauss-Seidel matrix solver
                phi(k) = GaussSeidel(A,b,Rtol);
            else
                phi(k) = A\b;
            end
        end
        for i = 1:N %Column sweep
            k = M*(i-1)+1:M*i;
            Yt = X; Yt(k,k) = 0;
            t = Yt*phi;
            b = S(k)-t(k);
            A = X(k,k);
            if strcmpi(method,'-TDMA') %use tridiagonal matrix solver
                phi(k) = TDMsolver(A,b);
            elseif strcmpi(method,'-GaussSeidel') %use Gauss-Seidel matrix solver
                phi(k) = GaussSeidel(A,b,Rtol);
            else
                phi(k) = A\b;
            end
        end        
        counter = counter + 1;
        R(counter+1) = norm(X*phi-S);
    end
    R = R(1:counter+1);
    phi = reshape(phi,M,N);
end