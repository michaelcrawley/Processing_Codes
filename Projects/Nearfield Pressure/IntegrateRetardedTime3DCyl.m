function p = IntegrateRetardedTime3DCyl(observer,field,source)
%field variable needs to be a struct with the following fields:
%       r   -   radial position     [M x N]
%       z   -   axial position      [M x N]
%       t   -   temporal position   [Q x 1]
%       c   -   average sound speed [1 x 1]
%source variable needs to be in the form [M x N x Q]

    %Constants
    Nphi = 4;
    phi = linspace(0,2*pi,Nphi+1);
    dphi = mean(diff(phi));

    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(field.z(1,:));
    [~,Ir] = sort(field.r(:,1));
    z = field.z(1,Iz);
    r = field.r(Ir,1);
    source = source(Ir,Iz,:);
    [M,N,T] = size(source);
    
    dt = abs(mean(diff(field.t)));
    dr = mean(diff(r));
    dz = mean(diff(z));
    
    voxel = repmat(pi*(r(2:end,1).^2 - r(1:end-1,1).^2),[1 N-1 Nphi])*dz*dphi/(2*pi); %integration volume
    voxel = voxel(:);
    source = repmat(permute(source,[1 2 4 3]),[1 1 Nphi+1 1]);    
    
    [Z,R,Phi] = meshgrid(z,r,phi);
    Z = Z(:);
    R = R(:);
    Phi = Phi(:);
    
    p = zeros(length(observer.z),T); 
%     pool = parpool;
    for k = 1:length(observer)
        D = sqrt((observer.r(k) - R.*cos(Phi)).^2 + (R.*sin(Phi)).^2 + (observer.z(k) - Z.^2)); %propagation distance
        tau = D/field.c; %propagation delay
        itau = round(tau/dt); %for now we are simply going to round to the nearest sampled time index, rather than interpolate
        D = reshape(D,[M N Nphi+1]);
%         Davg = (D(1:end-1,1:end-1,1:end-1) + D(2:end,1:end-1,1:end-1) + D(1:end-1,2:end,1:end-1) + D(1:end-1,1:end-1,2:end) + D(2:end,2:end,1:end-1) + ...
%                         D(2:end,1:end-1,2:end) + D(1:end-1,2:end,2:end) + D(2:end,2:end,2:end))/8;
        
        source = reshape(source,[],T);
        for Q = 1:M*N*(Nphi+1)
            source(Q,:) = circshift(source(Q,:),[0 itau(Q)]);
        end
        source = reshape(source,[M N Nphi+1 T]);
        
        for Q = 1:T
            source(:,:,:,Q) = source(:,:,:,Q)./D;
        end
        
        tstart = max(itau(:)); 
        for q = tstart:length(field.t)
            averaged = (source(1:end-1,1:end-1,1:end-1,q) + source(2:end,1:end-1,1:end-1,q) + source(1:end-1,2:end,1:end-1,q) + source(1:end-1,1:end-1,2:end,q) + source(2:end,2:end,1:end-1,q) + ...
                        source(2:end,1:end-1,2:end,q) + source(1:end-1,2:end,2:end,q) + source(2:end,2:end,2:end,q))/8;
                    
            integral = sum(averaged(:).*voxel);            
            p(k,q) = -integral/(4*pi*field.c*field.c);
        end
    end
%     delete(pool);
end