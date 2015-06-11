function M = perform_lic_v2(v, L, options)

% perform_lic - perform line integral convolution
%
%   M = perform_lic(v, L, options);
%
%   v is matrix (dimensions m by n by 2) which contains a vector field 
%     (should be approximately of unit norm).
%   L is the length of the convolution (in pixels)
%   M is an image of filtered noise that illustrates the flow of v.
%
%   options.spot_size is the desired width of blurs (in pixels)
%
%   Set options.flow_correction=1 in order to fix problem
%   around singularity points of the flow.
%
%   The method is described in
%       Imaging vector fields using line integral convolution
%       Brian Cabral, Leith Casey Leedom
%       Siggraph 1993
%  
%   
%   Copyright (c) Gabriel Peyre 2007
% Modified by Martin Kearney 2010 
%  -Converted to work with non-square vector fields.
%  -Removed broken code pieces.
%  -Reconfigured code for increased speed (roughly 30% faster) and memory handling.

n = size(v);
options.null = 0;   %Ensure struct "options" exists to prevent errors

    %Gets texture image if it exists, otherwise generates texture
if isfield(options, 'M0')   
    M0 = single(options.M0);
    options = rmfield(options,'M0');
else
    M0 = single(randn(n(1),n(2)));
    if isfield(options, 'spot_size')
        sigma = options.spot_size;
    else
        sigma = 2;
    end
    M0 = perform_blurring(M0, sigma, options);
end

    %Perform multiple iterations if requested - multiple iterations 
    %generally produce better results. Two is recommended for a good
    %balance between result quality and execution time.
if isfield(options, 'niter_lic') && options.niter_lic>1
    niter_lic = options.niter_lic;
else
    niter_lic = 1;
end

    %Gets or sets step-size for vector field integration (in pixels)
if isfield(options, 'dt')
    dt = options.dt;
else
    dt = 0.5;
end
    
    %Gets or sets histogram bins for histogram equalization
histogram = getoptions(options, 'histogram', 'gaussian');
if ischar(histogram)
    switch histogram
        case 'linear'
            hist = linspace(0,1, 100^2);
        case 'gaussian'
            hist = randn(100^2,1);
        otherwise
            disp('Unkown histogram type. Continuing without histogram equalization.');
            hist = [];
    end
else
    hist = histogram;
end   

verb = getoptions(options, 'verb', 0);  %Provide periodic progress updates

    %Perform integration of the vector field: Forward
T_list = [0 (0:dt:L/2)];
if verb~=0
    disp('Performing Forward Integration...');
end
H = perform_vf_integration(v, dt, T_list, options);

    %Perform integration of the vector field: Backward
T_list(1) = []; %Removes time equals zero to avoid double counting.
if verb~=0
    disp('Performing Backward Integration...');
end
H1 = perform_vf_integration(-v, dt, T_list, options); clear v;
H = cat(4, H1(:,:,:,end:-1:1), H); clear H1;
p = size(H,4);

    %Try to remove sampling problems
flow_correction = getoptions(options, 'flow_correction', 1);
if flow_correction
    A = H(:,:,:,(end+1)/2);
    dX = H(:,:,1,:)-repmat(A(:,:,1,:), [1 1 1 p]);
    dY = H(:,:,2,:)-repmat(A(:,:,2,:), [1 1 1 p]);  clear A;
    dX(dX>n(1)/2) = dX(dX>n(1)/2)-(n(1)-1); dX(dX<-n(1)/2) = dX(dX<-n(1)/2)+(n(1)-1);
    dY(dY>n(2)/2) = dY(dY>n(2)/2)-(n(2)-1); dY(dY<-n(2)/2) = dY(dY<-n(2)/2)+(n(2)-1);
    d = sqrt(dX.^2 + dY.^2); clear dX dY;
    d = repmat(mean(mean(d)), [n(1) n(2) 1] ) - d;

    eta = 0.5;  %Threshold on the distance map deviation
    
    VP = d < eta;   %Valid Points
    clear d;
else
    VP = ones(n(1),n(2),p,'single');
end


    %Compute averaging over streamlines for desired number of iterations
if verb~=0
    disp('Iterating...');
end
for j = 1:niter_lic  
    if verb~=0
        progress(j,1,niter_lic,verb);
    end
    M = zeros(n(1),n(2),'single');   N = M;  %N is a counter of the number of values added to a particular pixel
    for k = 1:p     %Iterate through each step in the streamline
        A = interp2(1:n(2),1:n(1), M0, H(:,:,2,k), H(:,:,1,k)); %Upscale texture image to match integrated result "H"
        M = M + A.*VP(:,:,k);
        N = N + VP(:,:,k);
    end
    N(N==0) = 1;    %Eliminates divide by zeros
    M = M./N;   %Computes average at each pixel

    if ~isempty(hist)
        M = perform_histogram_equalization(M, hist);
    end   
    
    M0 = M; %Sets texture image to current result for next iteration
end
 
