function [Ufilt,Vfilt,nmodes] = PODfilter(U,V,alpha)

    if ~exist('alpha','var')|| isempty(alpha), alpha = 0.999; end
    N = size(U);
    
    %Subtract out 'true' mean
    [~,Um] = nzstats(U,3);
    [~,Vm] = nzstats(V,3);
    Ufluct = U - repmat(Um,[1 1 N(3)]);
    Vfluct = V - repmat(Vm,[1 1 N(3)]);
    
    %Reshape and compute POD modes
    C = [reshape(Ufluct,prod(N(1:2)),N(3)); reshape(Vfluct,prod(N(1:2)),N(3))];
    [phi,lambda,ak] = SnapShotPOD(C,false);
    clear U V Ufluct Vfluct;
    
%     %Iteratively reconstruct filtered field, calculating the energy error
%     %as you go
%     deltaRM = zeros(N(3),1);
%     F = zeros(N(3)-2,1);
%     Uk = phi(:,1)*ak(1,:);
%     deltaRM(1) = std(Ufluct(:)-Uk(:));
%     for k = 2:N(3)
%         Ukp = phi(:,1:k)*ak(1:k,:);
%         deltaRM(k) = 
%     end

    %We're going to bullshit this....
    nmodes = round(N(3)/5);
    flag = true;
    while flag
        rec = phi(:,1:nmodes)*ak(1:nmodes,:);
        sigma = max(std(C-rec));
        F = (lambda(2:end)+N(3)*sigma^2)./(lambda(1:end-1)+N(3)*sigma^2);
        I = find(F > alpha,1) + 1;
        flag = I ~= nmodes;
        nmodes = I;
    end
    
    clear C;
    Ufilt = reshape(rec(1:end/2,:),N) + repmat(Um,[1 1 N(3)]);
    Vfilt = reshape(rec(end/2+1:end,:),N) + repmat(Vm,[1 1 N(3)]);
end