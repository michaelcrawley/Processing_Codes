function [F, N_c] = finite_diff_Taylor_uniform(N_pts,order,acc)
%   Assumptions:
%   1.  The field quantity whose derivative is desired is arranged as a
%       column vector with the successive elements storing the value at
%       successive grid points.
%   2.  The derivative would be obtained by premultiplying this column
%       vector by the square sparse matrix F returned from here.
%   3.  Since the operation is a matrix multiplication, multiple column
%       vectors of fields can be arranged in a 2D matrix for simultaneous
%       computation of the derivative.
%   4.  The grid is uniform.

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

if N_pts < min_grid_sz
    error(['Grid must have at least ', num2str(min_grid_sz), ' entries']);
end

%Compute the coefficients of the central difference scheme
%RHS of set of linear equations
rhs_c = zeros(K_c,1); 
rhs_c(order+1) = factorial(order); 
A_c = zeros(K_c);
for idx = 1:K_c
    A_c(idx,:) = (-N_c:N_c).^(idx-1);
end
coeffs_c = (A_c\rhs_c)';

%Temporarily create a sparse diagonal matrix with all points evaluated by
%the central-difference stencil.
F = spdiags(repmat(coeffs_c,N_pts,1), -N_c:N_c, N_pts, N_pts);

%Now fill in the coefficients for the forward and backward difference for
%the first and last N_c grid points (those that cannot be covered by the
%central difference). Note that these are anit--symmetric, so that only one
%end is computed.
%RHS of set of linear equations
rhs = zeros(K_a,1); 
rhs(order+1) = factorial(order); 
stencil = 1:K_a;
stencil_last = N_pts - (0:K_a-1);
for kdx = 1:N_c
    temp_stencil = stencil - kdx;
    A = zeros(K_a);
    for idx = 1:K_a
        A(idx,:) = temp_stencil.^(idx-1);
    end
    coeffs = A\rhs;
    F(kdx,stencil) = coeffs;
    F(N_pts-kdx+1,stencil_last) = -coeffs;
end
