function [F, N_c] = finite_diff_Taylor_non_uniform(r,order,acc)
%   Assumes that the field quantity (whose derivative is desired) is
%   arranged as a vector of the form
%       [q(r_1), ..., q(r_N)]'
%   where r_1...r_N are the coordinates of grid
%   
%   With the field quantity represented as above, the derivative operation
%   would be performed by simply pre-multiplying the corresponding vector
%   (or matrix) by the matrix F returned by this function.
%   Thus, the ith row of F represents the weighting of the appropriate
%   grid-points for obtaining the required order of finite difference for
%   the ith grid point, in the r-direction.
%   For interior points, the central difference is used. For boundary
%   points, the appropriate asymetric difference is used.
%   A complication that has to be handled here is the non-uniformity of the
%   r-grid.

%First determine the number of stencil points on either side of the pivot
%for central difference. Note that order of accuracy is even for central
%difference. So, the acc is increased to the closest higher even number if
%it is given as odd.
acc_c = 2*ceil(acc/2); %Accuracy in central difference
K_c = acc_c+2*ceil(order/2)-1; %Stencil size in central difference
N_c = round((K_c-1)/2); %Side points in central difference

%For all asymmetric finite differences, the stencil size is same. It is
%calculated here.
K_a = acc + order;

min_grid_sz = max(K_a, K_c);
N_r = length(r);

if N_r < min_grid_sz
    error(['r must have at least ', num2str(min_grid_sz), ' entries']);
end

F = zeros(N_r);
for kdx = 1:N_r
    if kdx <= N_c %Asymmetric forward difference
        stencil = 1:K_a;
    elseif kdx <= N_r - N_c %Central difference
        stencil = (-N_c:N_c)+kdx;
    else %Asymmetric backward difference
        stencil = N_r-K_a+1:N_r;
    end
    %RHS of set of linear equations
    rhs = zeros(length(stencil),1); 
    rhs(order+1) = factorial(order); 
    r_diff = r(stencil)-r(kdx);
    A = zeros(length(stencil));
    for idx = 1:length(stencil);
        A(idx,:) = r_diff.^(idx-1);
    end
    F(kdx,stencil) = A\rhs;
end
