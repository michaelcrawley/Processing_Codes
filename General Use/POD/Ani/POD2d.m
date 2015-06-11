function [phi, lambda] = POD2d(U)
%     L = length(flist);
% 
%     %read in first file to get sizing information, initialize variables
%     data = readimx([src_dir filesep flist{1}]);
%     [~,~,u,v] = showimx(data);
%     [M N] = size(u);
%     U = zeros(2*M*N,L); 
%     U(:,1) = [u(:); v(:)];    
%     
%     %read in data from all files
%     for n = 2:L
%         data = readimx([src_dir filesep flist{n}]);
%         [~,~,u,v] = showimx(data);        
%         U(:,n) = [u(:); v(:)];       
%     end
    [~,L] = size(U);
    
    %compute mean velocity vectors and subtract to get velocity
    %fluctuations
    umean = repmat(mean(U,2),1,L);
    U = U-umean;
        
    %compute inner product, find eigenvalues and eigenvectors
    corr = innerproduct(U,U)/L; %inner product is based on KE
    [~,eigenvalues,d] = svd(corr);
    lambda = diag(eigenvalues);
    
    %compute POD modes
    phi = U*d; 
    
    %scale individual POD modes
    for n = 1:L
        phi(:,n) = phi(:,n)/sqrt(innerproduct(phi(:,n),phi(:,n),c));
    end
    
    I = abs(lambda) <= 1000*eps;
    lambda(I) = 0;
    phi(:,I) = 0;
%     phi = reshape(phi,[size(u) L]);

end

function [R] = innerproduct(U,V)
    R = V'*U;
end