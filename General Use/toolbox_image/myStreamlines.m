function [X,Y,SL] = myStreamlines(xp,yp,u,v)
%This function creates an image of streamlines by blurring a texture image
%along the directions of the vector field. The function uses a third-party
%image processing toolbox called "Toolbox Image" and can be found at the
%following location (accessed: 2010-05-06):
%http://www.ceremade.dauphine.fr/~peyre/matlab/image/content.html

N = 4;    %Integer Upscaled Point density ratio - this is inversely related to the thickness of the streamlines.

% % % % u = CondAvg.U(xl(1):xl(2),yl(1):yl(2),1)-CondAvg.AvgUc; v = CondAvg.V(xl(1):xl(2),yl(1):yl(2),1);
% % % % xp = x(xl(1):xl(2),yl(1):yl(2)); yp = y(xl(1):xl(2),yl(1):yl(2));
%The LIC function requires a square field arranged as follows:
% 1) The u-component of the vector field should be in the first plane (e.g.
% V(:,:,1).
% 2) The v-component of the vector field should be in the second plane (e.g.
% V(:,:,2).
% 3) The complex functions used take a complex conjugate so the appropriate
% input complex function is: f(x,y) = u(x,y) +i*[-v(x,y)]. Therefore, v
% should be negated to produce proper streamlines.
% 4) For ease of computation, the vector field must be square. Therefore 
% the smaller dimension is padded with zeros.
S = size(u); n=max(S);
U = zeros(n);               %See above 4)
U(1:S(1),1:S(2)) = u;       %See above 1) & 4)
U(1:S(1),1:S(2),2) = -v;    %See above 2), 3), & 4)

    %The vector field is interpolated up to the desired resolution.
V = interp2(U(:,:,1),N); 
V(:,:,2) = interp2(U(:,:,2),N);  n = size(V,1); clear U;
V = perform_vf_normalization(V);    %The vector field is converted to unit vectors.

    % Input parameters for the LIC
LC = 25;                        %Length of convolution in pixels
options.bound = 'sym';          %Boundary handling - Is either 'per' (periodic extension) or 'sym' (symmetric extension).
options.histogram = 'linear';   %For contrast stretching - Is either 'linear' (keep contrast fixed - just stretch) or 'gaussian' (use gaussian contrast distribution).
options.verb = 0;               %Logical to show progress bar (1==show, 0==don't show).
options.spot_size = 1;          %The width of blurs (in pixels)
options.dt = 1.5;               %Step-size for vector field integration (in pixels)
options.flow_correction = 1;    %Correcting the flow (==1) removes singularities
options.niter_lic = 2;          %Iterate N times through the convolution process. Multiple iterations generally give better results

if n > 800
    OL = 0.05;  %Overlap percentage
    ns = ceil(n/((1-OL)*800));      %Number of segments in one dimension
    nps = ceil(n/(ns*(1-OL)+OL));   %Number of points per segment
    nOL = ceil(OL*nps); %Number of points overlap
    
    np = nps + (ns-1)*(nps-nOL);    %Number of points 
    V(np,np,:) = 0;
    
    IC = [1 (nps-nOL+1:nps-nOL:np-nOL-1)];
    IC(2,:) = IC+nps-1;
    
    TI = randn(nps);   %Creates the texture image of white noise - white noise generally works best.
    options.M0 = TI; clear TI;      %Add texture image to processing options
    
    VM = zeros(nps); SL = zeros(np); SM = SL; I = 0;
    for j = 1:ns
        for k = 1:ns
            I = I+1;
            fprintf(['..Image Segment ' num2str(I) ' of ' num2str(ns^2) '\n']);
            VM = V(IC(1,j):IC(2,j),IC(1,k):IC(2,k),:);
            SL(IC(1,j):IC(2,j),IC(1,k):IC(2,k)) = SL(IC(1,j):IC(2,j),IC(1,k):IC(2,k)) +perform_lic(VM, LC, options);   %Compute streamlines
            SM(IC(1,j):IC(2,j),IC(1,k):IC(2,k)) = SM(IC(1,j):IC(2,j),IC(1,k):IC(2,k)) +1;
        end
    end
    SL = SL./SM;   
        
else
    TI = randn(n);   %Creates the texture image of white noise - white noise generally works best.
    options.M0 = TI; clear TI;      %Add texture image to processing options

    SL = perform_lic(V, LC, options);   %Compute streamlines
end

clear V options
X = interp2(xp,N); Y = interp2(yp,N); Sp = size(X); %Upscale coordinates for streamline image
SL = SL(1:Sp(1),1:Sp(2)); clear S Sp;   %Crop off zero padding

%% PLOTTING %%%%%%%%%%%%%%%%%%%%%
% Q = CondAvg.L(xl(1):xl(2),yl(1):yl(2),1);
% L = interp2(Q/max(Q(:)),N);
% NC = 5; CM = makeLine([0 1 1],[1 0 1],NC);
% L = round(L*size(CM,1));
% I = L==0;
% Lrgb = ind2rgb(L,CM);
% Lycbcr = rgb2ycbcr(Lrgb);
% Lycbcr(:,:,1) = SL;
% Lrgb2 = ycbcr2rgb(Lycbcr);
% tmp = Lrgb2(:,:,1); tmp(I) = M(I); Lrgb2(:,:,1) = tmp;
% tmp = Lrgb2(:,:,2); tmp(I) = M(I); Lrgb2(:,:,2) = tmp;
% tmp = Lrgb2(:,:,3); tmp(I) = M(I); Lrgb2(:,:,3) = tmp;
% 
% P = permute(Lrgb2,[2 1 3]);
% image(X(:,1),Y(1,:),P)
% set(gca,'YDir','normal')
% axis([1 8 -2 2])
% colormap(CM);
% colorbar;
% CT = round((max(Q(:))-min(Q(:)))/NC*((1:NC)-0.5));  YTL = cell(1,NC);
% for n = 1:NC
%     YTL{n} = num2str(CT(n));
% end
% colorbar('YTick',(1:NC),'YTickLabel',YTL)
% xlabel('x/D')
% ylabel('y/D')

