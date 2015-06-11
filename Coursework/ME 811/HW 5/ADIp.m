function [phi,R] = ADIp(X,S,method,Rtol,itrmax)
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

    AW = full([zeros(M,1); diag(X,-M)]);
    AE = full([diag(X,M); zeros(M,1)]);
    AN = full([diag(X,1); 0]);
    AS = full([0; diag(X,-1)]);
    AO = full(diag(X));
    
    while (R(counter+1) >= Rtol) && (counter < itrmax)        
        %Row sweep
        k = 1:M:1+M*(N-1);
        b = S(k)-AN(k).*phi(k+1);
        phi(k) = TDMAsolver(AW(k),AO(k),AE(k),b);
        for j = 2:M-1
            k = j:M:j+M*(N-1); %determine nodal array points
            b = S(k)-AN(k).*phi(k+1)-AS(k).*phi(k-1);
            phi(k) = TDMAsolver(AW(k),AO(k),AE(k),b);  
        end
        k = M:M:M+M*(N-1); %determine nodal array points
        b = S(k)-AS(k).*phi(k-1);
        phi(k) = TDMAsolver(AW(k),AO(k),AE(k),b);
        
        %Column sweep
        k = M*(1-1)+1:M*1;        
        b = S(k)-AE(k).*phi(k+M);
        phi(k) = TDMAsolver(AS(k),AO(k),AN(k),b);
        for i = 2:N-1 
            k = M*(i-1)+1:M*i;
            b = S(k)-AE(k).*phi(k+M)-AW(k).*phi(k-M);
            phi(k) = TDMAsolver(AS(k),AO(k),AN(k),b);
        end 
        k = M*(N-1)+1:M*N;
        b = S(k)-AW(k).*phi(k-M);
        phi(k) = TDMAsolver(AS(k),AO(k),AN(k),b);        
        
        counter = counter + 1;
        R(counter+1) = norm(X*phi-S);
    end

    R = R(1:counter+1);
    phi = reshape(phi,M,N);
end