function [X,Y,SL] = myStreamlines(x,y,u,v)
%This function creates an image of streamlines by blurring a texture image
%along the directions of the vector field. The function uses a third-party
%image processing toolbox called "Toolbox Image" and can be found at the
%following location (accessed: 2010-05-06):
%http://www.ceremade.dauphine.fr/~peyre/matlab/image/content.html

%The fast code implementation is about 30% faster than the low memory code,
%but its memory consumption is heavily dependent upon: "size(x)", "N", 
%"LC", and "dt". The memory use of the low memory consumption code depends
%only on "size(x)" and "N".

    % Input parameters for the LIC
N = 3;    %Integer Upscaled Point density ratio - this is inversely related to the thickness of the streamlines.
%   4 is usually a good choice
FL = 'fast';    %Fast code or low memory consumption code - 'fast' or 'low'
LC = 30;                        %Length of convolution in pixels
options.bound = 'sym';          %Boundary handling - Is either 'per' (periodic extension) or 'sym' (symmetric extension).
options.histogram = 'linear';   %For contrast stretching - Is either 'linear' (keep contrast fixed - just stretch) or 'gaussian' (use gaussian contrast distribution).
options.verb = 10;               %Logical to show progress updates (0==don't show, else==update interval).
options.dt = 1;               %Step-size for vector field integration (in pixels)
options.flow_correction = 1;    %Correcting the flow (==1) removes singularities
options.spot_size = 3;          %The width of blurs (in pixels)
options.niter_lic = 2;          %Iterate N times through the convolution process. Multiple iterations generally give better results
% options.M0 = TI; clear TI;      %Add your own texture image to processing options - leaving "options.M0" non-existent causes code to use a random texture with spot sizes specified by "options.spot_size"

    %Rotate data into proper orientation
Fflags = [0 0 0];   %Matrix manipulation flags: Transpose, FlipUD, FlipLR
if diff(x(1:2,1))~=0    
    x = x'; y = y'; u = u'; v = v';
    Fflags(1) = 1;
end
if diff(y(1:2,1))<0
    x = flipud(x); y = flipud(y); u = flipud(u); v = flipud(v);
    Fflags(2) = 1;
end
if diff(x(1,1:2))<0
    x = fliplr(x); y = fliplr(y); u = fliplr(u); v = fliplr(v);
    Fflags(3) = 1;
end    

%The LIC function requires a square field arranged as follows:
% 1) The u-component of the vector field should be in the second plane (e.g.
% V(:,:,2).
% 2) The v-component of the vector field should be in the first plane (e.g.
% V(:,:,1).
% 3) The complex functions use the following form for the streamfunction: 
%        f(x,y) = u(x,y) +i*v(x,y).

    %The vector field is interpolated up to the desired resolution.
V = interp2(single(v),N); 
V(:,:,2) = interp2(single(u),N);
V = normalize_vf(V);    %The vector field is converted to unit vectors.

if strcmpi(FL,'fast')   %Compute streamlines
    disp('Using Fast Code...')
    SL = perform_lic_v2(V, LC, options);
else
    disp('Using Low Memory Code...')
    SL = perform_lic_v3(V, LC, options);   
end    
clear V options;

    %Return outputs to original orientation
if Fflags(3)
    x = fliplr(x); y = fliplr(y); SL = fliplr(SL);
end
if Fflags(2)
    x = flipud(x); y = flipud(y); SL = flipud(SL);
end
if Fflags(1)
    x = x'; y = y'; SL = SL';
end
X = interp2(x,N); Y = interp2(y,N); %Upscale coordinates for streamline image
