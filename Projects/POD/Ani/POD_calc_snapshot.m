function [S, Phi, U_mean, C] = POD_calc_snapshot(u, inner_product, inner_product_opts, flag_mean, N_modes)

% Calculates POD modes by snapshot method.
% Parameters:
%     Input:
%         u                   : The field. M columns for M snapshots. The
%                               rows would typically be arranged as
%                               follows for 3D velocity defined on a 2d
%                               NX-by-NY grid:
%                               U(:,t) = [  u(1,1,t) ... u(NX,1,t) ... u(1,NY,t) ... u(NX,NY,t) 
%                                           v(1,1,t) ... v(NX,1,t) ... v(1,NY,t) ... v(NX,NY,t)
%                                           w(1,1,t) ... w(NX,1,t) ... w(1,NY,t) ... w(NX,NY,t) ]'
%         inner_product       : Inner product function reference
%         inner_product_opts  : Parameters needed by above
%         flag_mean           : 1 = Subtract mean, 0 = Do not.
%         N_modes             : Number of POD modes desired
%     Output:
%         C         :   Cross-correlation of snapshots 
%         S         :   All M eigenvalues
%         Phi       :   POD modes. N rows for N grid points, N_modes
%                       columns for N_modes POD modes
%         U_mean    :   Mean field for flag_mean = 1, 0-vector otherwise.

[N,M] = size(u);

switch nargin
    case 1
        inner_product = @(a,b,c)(b'*diag(c)*a);
        inner_product_opts = ones(N,1);
        flag_mean = 0;
        N_modes = M;
    case 2
        inner_product_opts = ones(N,1);
        flag_mean = 0;
        N_modes = M;
    case 3
        flag_mean = 0;
        N_modes = M;
    case 4
        N_modes = M;        
    case 5
    otherwise
        error('Incorrect number of arguments')
end

if isempty(inner_product)
    inner_product = @(a,b,c)(b'*diag(c)*a);
end
if isempty(inner_product_opts)
    inner_product_opts = ones(N,1);
end
if isempty(flag_mean)
    flag_mean = 0;
end

if flag_mean
    U_mean = mean(u,2);
    u = u - repmat(U_mean,1,M);
else
    U_mean = zeros(N,1);
end

C = inner_product(u,u,inner_product_opts)/M;

[U,S,V] = svd(C);
S = diag(S);

Phi = zeros(N,N_modes);
for idx = 1:N_modes
    Phi(:,idx) = u*V(:,idx);    
    Phi(:,idx) = Phi(:,idx)/sqrt(inner_product(Phi(:,idx),Phi(:,idx),inner_product_opts));  %Normalization
end
