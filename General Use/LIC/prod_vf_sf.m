function v2 = prod_vf_sf(v1,s)
%Compute the product of a vector field by a scalar field.

v2 = zeros(size(v1));

for i=1:size(v1,3)
    v2(:,:,i) = v1(:,:,i).*s;
end