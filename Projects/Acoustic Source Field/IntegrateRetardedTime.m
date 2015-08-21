function p = IntegrateRetardedTime(observer,field,source)
%field variable needs to be a struct with the following fields:
%       r   -   radial position     [M x N]
%       z   -   axial position      [M x N]
%       t   -   temporal position   [Q x 1]
%       c   -   average sound speed [1 x 1]
%source variable needs to be in the form [M x N x Q]

    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(field.z(1,:));
    [~,Ir] = sort(field.r(:,1));
    z = field.z(1,Iz);
    r = field.r(Ir,1);
    source = source(Ir,Iz,:);
    [M,N,~] = size(source);
    
    dt = abs(mean(diff(field.t)));
    dr = mean(diff(r));
    dz = mean(diff(z));
     
    p = cell(length(observer.z),1); 

    for k = 1:length(observer)
        R = sqrt((observer.r(k) - field.r).^2 + (observer.z(k) - field.z).^2); %propagation distance
        tau = R./field.c; %propagation delay
        itau = round(tau/dt); %for now we are simply going to round to the nearest sampled time index, rather than interpolate
        
        tstart = max(itau(:));
        
        retarded_source = zeros(size(source));
        for m = 1:M
            for n = 1:N
                retarded_source(m,n,:) = circshift(source(m,n,:),[0 0 itau(m,n)]);
            end
        end
        clear source;
        
        
        for q = tstart:length(field.t)
            integral = trapz(trapz(retarded_source(:,:,q),1)*dr,2)*dz;
            p{k}(q) = -integral/(4*pi*field.c*field.c);
        end
    end
end