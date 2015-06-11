function [x,R] = GaussSeidel(A,b,Rtol,itrmax)
%   Solves linear system of equations using the Gauss-Seidel method.
%Inputs:
%       A:      Full matrix, where Ax=b
%       b:      Column vector
%       Rtol:   Residual tolerance (optional; default = 10E-6)
%       itrmax: Maximum number of iterations for the algorithm to perform
%               (optional; default = 10,000)
%Outputs:
%       x:      Solution to the linear equation
%       R:      Final residual


    if ~exist('Rtol','var')
        Rtol = 10E-6;
    end
    if ~exist('itrmax','var')
        itrmax = 10000;
    end
    x = zeros(size(b));
   
    R = sqrt(sum((A*x-b).^2));
    itr = 1;
    while R > Rtol && itr <= itrmax
        for i = 1:length(x)
           x(i) = (1/A(i,i))*(b(i)-A(i,i+1:end)*x(i+1:end)-A(i,1:i-1)*x(1:i-1)); 
        end
        R = sqrt(sum((A*x-b).^2));
        itr = itr + 1;
        if any(isnan(x)) || any(isinf(x))
            warning('Error: Gauss-Seidel Method failed to converge'); %#ok<WNTAG>
            break;
        end        
    end
    
    if itr == itrmax && R > Rtol 
        warning(['Error: Gauss-Seidel Method failed to converge within ',num2str(itrmax),' iterations']); %#ok<WNTAG>
    end
end