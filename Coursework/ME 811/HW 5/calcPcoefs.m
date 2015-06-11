function [ApO ApE ApW ApN ApS Sp] = calcPcoefs(alpha,rho,dx,dy,uh,vh,AxO,AxE,AxW,AxN,AxS,AyO,AyE,AyW,AyN,AyS)
    
    [M,~] = size(AxO);
    [~,N] = size(AyO);
    
    sAx = [Inf*ones(M,1) ((1+alpha)*AxO+AxE+AxW+AxN+AxS) Inf*ones(M,1)]; %AxO = (:,1:end-1), AxE = (:,2:end) 
    sAy = [Inf*ones(1,N);  ((1+alpha)*AyO+AyE+AyW+AyN+AyS); Inf*ones(1,N)]; %AyO = (2:end,:), AyN = (1:end-1,:)
    uh = [uh zeros(M,1)];
    vh = [zeros(1,N);  vh];
    
    %Calculate link coefficients for all nodes
    ApE = rho*dy^2./sAx(:,2:end);
    ApW = rho*dy^2./sAx(:,1:end-1);
    ApN = rho*dx^2./sAy(1:end-1,:);
    ApS = rho*dx^2./sAy(2:end,:);
    ApO = -(ApE+ApW+ApN+ApS);
    Sp = rho*((uh(:,2:end)-uh(:,1:end-1))*dy+(vh(1:end-1,:)-vh(2:end,:))*dx);
end