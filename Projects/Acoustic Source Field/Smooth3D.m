function smoothed = Smooth3D(signal,whan,stencil)

    [M,N,K] = size(signal);
    num = numel(signal);
    L = (stencil-1)/2;
    DIM1_coefs = whan(stencil(1))'/sum(whan(stencil(1)));
    DIM2_coefs = whan(stencil(2))'/sum(whan(stencil(2)));
    DIM3_coefs = whan(stencil(3))'/sum(whan(stencil(3)));
    
    sparse = spdiags(repmat(DIM1_coefs,num,1),-L:L,num,num) + spdiags(repmat(DIM2_coefs,num,1),-(L*M):M:(L*M),num,num) + spdiags(repmat(DIM3_coefs,num,1),-(L*M*N):(M*N):(L*M*N),num,num);
    signal = reshape(signal,[],1);
    smoothed = sparse*signal;
    smoothed = reshape(smoothed,[M,N,K]);
end