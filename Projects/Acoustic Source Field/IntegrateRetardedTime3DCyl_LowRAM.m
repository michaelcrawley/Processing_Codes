function p = IntegrateRetardedTime3DCyl_LowRAM(observer,field,source)
%field variable needs to be a struct with the following fields:
%       r   -   radial position     [M x N]
%       z   -   axial position      [M x N]
%       t   -   temporal position   [Q x 1]
%       c   -   average sound speed [1 x 1]
%source variable needs to be in the form [M x N x Q]

    %Constants
    Nphi = 8;
    phi = linspace(0,2*pi,Nphi+1);
    dphi = mean(diff(phi));
    
    %Make sure matrices are in correct ascending order
    [~,Iz] = sort(field.z(1,:));
    [~,Ir] = sort(field.r(:,1));
    z = field.z(1,Iz);
    r = field.r(Ir,1);
    source = source(Ir,Iz,:);
    [M,N,T] = size(source);
    
    pool = parpool(12);
    %Set up progress bar
    cpb = ConsoleProgressBar();
    cpb.setLeftMargin(4);   % progress bar left margin
    cpb.setTopMargin(1);    % rows margin
    cpb.setLength(40);      % progress bar length: [.....]
    cpb.setMinimum(0);      % minimum value of progress range [min max]
    cpb.setMaximum(100);    % maximum value of progress range [min max]
    cpb.setPercentPosition('left');
    cpb.setTextPosition('right');
    cpb.setMinimum(0);
    cpb.setMaximum(T);
    cpb.start();
    
    dt = abs(mean(diff(field.t)));
    dz = mean(diff(z));
    
    voxel = repmat(pi*(r(2:end,1).^2 - r(1:end-1,1).^2),[1 N-1 Nphi])*dz*dphi/(2*pi); %integration volume
    voxel = voxel(:);  
    
    [Z,R,Phi] = meshgrid(z,r,phi);
    
    p = zeros(length(observer.z),T); 
    
    for k = 1:length(observer)
        D = sqrt((observer.r(k) - R.*cos(Phi)).^2 + (R.*sin(Phi)).^2 + (observer.z(k) - Z.^2)); %propagation distance
        tau = D/field.c; %propagation delay
        itau = round(tau/dt); %for now we are simply going to round to the nearest sampled time index, rather than interpolate              
                
        tstart = max(itau(:))+1;         
        for q = tstart:length(field.t)
            retarded_time = zeros(M*N*(Nphi+1),1);
            local = source(:,:,q-tstart+1:q);
            parfor Q = 1:M*N*(Nphi+1)
                [m,n,qq] = ind2sub([M N (Nphi+1)],Q);
                retarded_time(Q) = local(m,n,end-itau(m,n,qq));
            end
            retarded_time = reshape(retarded_time,[M,N,Nphi+1]);
            
            averaged = (retarded_time(1:end-1,1:end-1,1:end-1) + retarded_time(2:end,1:end-1,1:end-1) + retarded_time(1:end-1,2:end,1:end-1) + retarded_time(1:end-1,1:end-1,2:end) + retarded_time(2:end,2:end,1:end-1) + ...
                        retarded_time(2:end,1:end-1,2:end) + retarded_time(1:end-1,2:end,2:end) + retarded_time(2:end,2:end,2:end))/8;
                    
            integral = voxel'*averaged(:);            
            p(k,q) = -integral/(4*pi*field.c*field.c);
             
            text = sprintf('Iteration: %d/%d', q, T);
            cpb.setValue(q);  	% update progress value
            cpb.setText(text);  % update user text
        end
    end
    delete(pool);
end