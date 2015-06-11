function v2 = normalize_vf(v1)
%Normalize a vector field.


n = sum(v1.^2,3);
I = n<eps;  %Prevents NaN creation due to divide by zero.
n(I) = 1;
v2 = prod_vf_sf(v1,1./sqrt(n));