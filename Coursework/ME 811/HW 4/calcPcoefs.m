function [ApO ApE ApW ApN ApS Sp] = calcPcoefs(rho,dx,dy,uh,vh,AxO,AyO)
    
    [M,~] = size(AxO);
    [~,N] = size(AyO);
    
    AxO = [Inf*ones(M,1) AxO Inf*ones(M,1)]; %AxO = (:,1:end-1), AxE = (:,2:end) 
    AyO = [Inf*ones(1,N);  AyO; Inf*ones(1,N)]; %AyO = (2:end,:), AyN = (1:end-1,:)
    uh = [uh zeros(M,1)];
    vh = [zeros(1,N);  vh];
    
    %Calculate link coefficients for all nodes
    ApO = -rho*(dy^2./AxO(:,2:end)+dy^2./AxO(:,1:end-1)+dx^2./AyO(1:end-1,:)+dx^2./AyO(2:end,:));
    ApE = rho*dy^2./AxO(:,2:end);
    ApW = rho*dy^2./AxO(:,1:end-1);
    ApN = rho*dx^2./AyO(1:end-1,:);
    ApS = rho*dx^2./AyO(2:end,:);
    Sp = rho*((uh(:,2:end)-uh(:,1:end-1))*dy+(vh(1:end-1,:)-vh(2:end,:))*dx);
end