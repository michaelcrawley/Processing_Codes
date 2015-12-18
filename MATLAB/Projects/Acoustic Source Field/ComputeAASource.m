function S = ComputeAASource(r,z,rho,Ur,Uz)
    [M,N] = size(Uz);
    %First DIM should be r, second DIM z
    
    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(z(1,:));
    [~,Ir] = sort(r(:,1));
    z = z(:,Iz);
    r = r(Ir,:);
    Uz = Uz(Ir,Iz);
    Ur = Ur(Ir,Iz);
        
    %Create derivative matrices for divergence calculations
    order = 4;
    dz = mean(diff(z(1,:)));
    dr = mean(diff(r(:,1))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    partial_z = mNumericalDerivative(1,order,dz,N);
    partial_r = mNumericalDerivative(1,order,dr,M);
    
    %Compute uXu tensor
    Tensor{1} = rho.*Ur.*Ur; %T_rr
    Tensor{2} = rho.*Ur.*Uz; %T_rz = T_zr
    Tensor{3} = rho.*Uz.*Uz; %T_zz 
    
    %Compute divergence of tensor (vector)
    Vector{1} = partial_r*Tensor{1} + Tensor{1}./r + (partial_z*Tensor{2}')';
    Vector{2} = partial_r*Tensor{2} + Tensor{2}./r + (partial_z*Tensor{3}')';
    
    %Compute Source field
    S = partial_r*Vector{1} + Vector{1}./r + (partial_z*Vector{2}')';
end