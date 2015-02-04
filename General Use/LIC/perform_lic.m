function M = perform_lic(v, w, options)

% perform_lic - perform line integral convolution
%
%   M = perform_lic(v, w, options);
%
%   v is a vector field (should be approximately of unit norm).
%   w is the length of the convolution (in pixels)
%   M is an image of filtered noise that illustrate well the flow of v.
%
%   options.spot_size set the size of the features.
%
%   Set options.flow_correction=1 in order to fix problem
%   around singularity points of the flow.
%
%   The method is described in
%       Imaging vector fields using line integral convolution
%       Brian Cabral, Leith Casey Leedom
%       Siggraph 1993
%
%   If the vector field is not oriented (e.g. if it is an eigenvector field
%   of a tensor field) then set options.isoriented=1. It will perform two
%   LIC (horizontal and vertical) and then average them. This trick is
%   explained in 
%       Interactive Tensor Field Design and Visualization on Surfaces 
%       Eugene Zhang, James Hays, and Greg Turk 
%       IEEE Transactions on Visualization and Computer Graphics
%       Volume 13 ,  Issue 1, Pages: 94-107, 2007.   
%   
%   Copyright (c) Gabriel Peyre 2007
% Modified by Martin Kearney 2010 
%  -Converted to work with non-square vector fields.
%  -Removed broken code pieces.

n = size(v);

options.null = 0;
if isfield(options, 'M0')
    M0 = options.M0;
else
    M0 = randn(n(1),n(2));
    if isfield(options, 'spot_size')
        sigma = options.spot_size;
    else
        sigma = 2;
    end
    M0 = perform_blurring(M0, sigma, options);
end

    %Perform multiple iterations if requested
if isfield(options, 'niter_lic') && options.niter_lic>1
    M = options.M0;
    niter_lic = options.niter_lic;
    options.niter_lic = 1;
    for i=1:niter_lic
        options.M0 = M;
        M = perform_lic(v, w, options);
    end
    return;
end


histogram = getoptions(options, 'histogram', 'gaussian');
flow_correction = getoptions(options, 'flow_correction', 1);

if ischar(histogram)
    switch histogram
        case 'linear'
            hist = linspace(0,1, 100^2);
        case 'gaussian'
            hist = randn(100); hist = hist(:);
        otherwise
            error('Unkown kind of histograms');
    end
else
    hist = histogram;
end

if isfield(options, 'dt')
    dt = options.dt;
else
    dt = 0.5;
end

% perform integration of the vector field: Forward
T_list = 0:dt:w/2;
H = perform_vf_integration(v, dt, T_list, options );

% perform integration of the vector field: Backward
T_list(1) = [];
H1 = perform_vf_integration(-v, dt, T_list, options );
H = cat(4, H1(:,:,:,end:-1:1), H ); clear H1;
p = size(H,4);

% try to remove sampling problems
if flow_correction
    A = H(:,:,:,(end+1)/2);
    dX = H(:,:,1,:)-repmat(A(:,:,1,:), [1 1 1 p]);
    dY = H(:,:,2,:)-repmat(A(:,:,2,:), [1 1 1 p]);
    dX(dX>n(1)/2) = dX(dX>n(1)/2)-(n(1)-1); dX(dX<-n(1)/2) = dX(dX<-n(1)/2)+(n(1)-1);
    dY(dY>n(2)/2) = dY(dY>n(2)/2)-(n(2)-1); dY(dY<-n(2)/2) = dY(dY<-n(2)/2)+(n(2)-1);
    d = sqrt(dX.^2 + dY.^2);
    d = repmat( mean(mean(d)), [n(1) n(2) 1] ) - d;
    
    eta = 0.5;  % threshold on the distance map deviation
end

% compute averaging over streamline
M = zeros(n(1),n(2));   W = M;
for i=1:p
    A = interp2(1:n(2),1:n(1), M0, H(:,:,2,i), H(:,:,1,i) );
    if flow_correction
        w = d(:,:,i)<eta;
    else
        w = ones(n(1),n(2));
    end
    M = M + A.*w;
    W = W+w;
end
W(W==0) = 1;
M = M./W;

if not(isempty(hist))
    M = perform_histogram_equalization(M, hist);
end
    
