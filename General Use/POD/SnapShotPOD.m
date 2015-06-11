function [phi,lambda,ak] = SnapShotPOD(U,flag)

    %Using notation set forth in Meyer, "Proper othrogonal decomposition on
    %a jet in crossflow"
    %DIM1 of U is position, DIM2 of U is snapshot
    
    if ~exist('flag','var')||isempty(flag), flag=true; end
    
    N = size(U);
    
    if flag
        U = U - repmat(mean(U,2),1,N(2)); %subtract mean value at each location
    end
    C = U.'*U; %correlation matrix in space
    [~,eigV,A] = svd(C);
    lambda = diag(eigV);
    phi = U*A;
    normal = sqrt(sum(phi.^2,1)); %sum over spatial locations
    phi = phi./repmat(normal,N(1),1); %normalize POD modes
    ak = phi.'*U;
end