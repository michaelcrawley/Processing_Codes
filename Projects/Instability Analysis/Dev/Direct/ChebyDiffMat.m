function [D] = ChebyDiffMat(N,ndiff)
%Creates differentiation matrices for use in Chebyshev collocation method.
%Inputs:
%       N:      Number of Chebyshev modes
%       ndiff:  Highest derivative order required
%
%Outputs:
%       D:      (N+1,N+1,ndiff+1) matrix containing all the Chebyshev
%       differentiation matrices (order is (:,:,1:ndiff+1)).

    D = zeros(N+1,N+1,ndiff+1);
    D(:,:,1) = cos(((0:N)'*(0:N))*(pi/N));
    
    if ndiff >= 2
        D(:,1:3,2) = [zeros(N+1,1) D(:,1,1) 4*D(:,2,1)];
        D(:,1:3,3) = [zeros(N+1,2) 4*D(:,1,1)];
    elseif ndiff >= 1
        D(:,1:3,2) = [zeros(N+1,1) D(:,1,1) 4*D(:,2,1)];
    end
    
    for j = 3:N
       D(:,j+1,2:end) = 2*j*D(:,j,1:end-1)+j*D(:,j-1,2:end)/(j-2);
    end
end