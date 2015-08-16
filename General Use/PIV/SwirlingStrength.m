function L = SwirlingStrength(varargin)
%This program calculates the \Delta-Criterion (D), the Swirling Strength
%(L), the Q-Criterion (Q), and the \lambda_2-Criterion (L2) based on an
%input flow field. These methods require incompressibility to work
%properly. In the event that only two velocity components are provided, the
%program enforces the incompressibility condition. To clarify,
%mathematically speaking, incompressibility is required, but, in reality,
%quasi-incompressible is acceptable.
%
%Valid calls:
% [D,L,Q,L2] = VortexID(u,v);
%    -Grid step-size set to 1 and w-component computed to satisfy
%    incompressibility.
% [D,L,Q,L2] = VortexID(u,v,w);
%    -Grid step-size set to 1.
% [D,L,Q,L2] = VortexID(x,y,u,v);
%    -w-component computed to satisfy incompressibility.
% [D,L,Q,L2] = VortexID(x,y,z,u,v,w);
%
%INPUTS
% x - x-axis coordinates - should be of same dimensions as velocity field
%   such as would be formed by meshgrid.
% y - y-axis coordinates - same requirement as "x"
% z - z-axis coordinates - same requirement as "x"
% u - x-velocity component
% v - y-velocity component
% w - z-velocity component
%
%OUTPUTS
% D - \Delta-Criterion
% L - Swirling Strength (1\sec)
% Q - Q-Criterion
% L2 - \lambda_2-Criterion
%
%Written by: Martin Kearney-Fischer - 05-03-2010

%% INPUT Formatting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag2D = false;
dlo = [2 1 3];  %matrix dimension along which x, y, and z vary respectively for calculations
switch nargin
    case 2    %VortexID(u,v)
        flag2D = true;
    
        u = varargin{1};
        v = varargin{2};

            %Creates coordinates and missing velocity component so that
            %derivatives may be easily calculated.
        u = repmat(u,[1 1 3]);
        v = repmat(v,[1 1 3]);
        [x,y,z] = meshgrid((1:size(u,1)),(1:size(u,2)),[0 1 2]);
        w = zeros(size(u));

        dli = dlo;  %matrix dimension along which x, y, and z vary respectively
    case 3    %VortexID(u,v,w)
        u = varargin{1};
        v = varargin{2};
        w = varargin{3};
        
        [x,y,z] = meshgrid((1:size(u,1)),(1:size(u,2)),(1:size(u,3)));

        dli = dlo;  %matrix dimension along which x, y, and z vary respectively
    case 4    %VortexID(x,y,u,v)
        flag2D = true;

        x = varargin{1};
        y = varargin{2};
        u = varargin{3};
        v = varargin{4};

        if x(2,1)-x(1,1)==0
            [x,y,z] = meshgrid(x(1,:),y(:,1),[0 1 2]);
            u = repmat(u,[1 1 3]);
            v = repmat(v,[1 1 3]);

            dli = dlo;  %matrix dimension along which x, y, and z vary respectively
        else
            [x,y,z] = meshgrid(x(:,1),y(1,:),[0 1 2]);
            u = u'; u = repmat(u,[1 1 3]);
            v = v'; v = repmat(v,[1 1 3]);

            dli = [1 2 3];  %matrix dimension along which x, y, and z vary respectively
        end
        w = zeros(size(u));
    case 6  %VortexID(x,y,z,u,v,w)
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        u = varargin{4};
        v = varargin{5};
        w = varargin{6};

            %Determines variation direction for coordinates
        if (x(2,1,1)-x(1,1,1)==0) && (x(1,2,1)-x(1,1,1)==0)
            dli(1) = 3;
        elseif (x(1,1,2)-x(1,1,1)==0) && (x(1,2,1)-x(1,1,1)==0)
            dli(1) = 1;
        else
            dli(1) = 2;
        end
        if (y(2,1,1)-y(1,1,1)==0) && (y(1,2,1)-y(1,1,1)==0)
            dli(2) = 3;
        elseif (y(1,1,2)-y(1,1,1)==0) && (y(1,2,1)-y(1,1,1)==0)
            dli(2) = 1;
        else
            dli(2) = 2;
        end
        if (z(2,1,1)-z(1,1,1)==0) && (z(1,2,1)-z(1,1,1)==0)
            dli(3) = 3;
        elseif (z(1,1,2)-z(1,1,1)==0) && (z(1,2,1)-z(1,1,1)==0)
            dli(3) = 1;
        else
            dli(3) = 2;
        end
        if length(unique(dli))~=3
            error('Poorly formed coordinates')
        end

        if sum(dli==dlo)~=3 %Reorders the matrices to conform to calculation method
            cx = zeros(3,1);
            for n = 1:3
                cx(n) = find(dli==dlo(n));
            end
            
            x = permute(x,cx);
            y = permute(y,cx);
            z = permute(z,cx);
            u = permute(u,cx);
            v = permute(v,cx);
            w = permute(w,cx);
        end
    otherwise
    error('Incorrect number of input arguments');
end

%% CORE PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = size(u); L = zeros(S);     %Initializes output arrays

    %Calculates second order accurate derivatives for all components
dudx = zeros(size(u)); dvdx = dudx; dwdx = dudx; dudy = dudx; 
dvdy = dudx; dwdy = dudx; dudz = dudx; dvdz = dudx; dwdz = dudx;
dudx(:,2:end-1,:) = (u(:,3:end,:)-u(:,1:end-2,:))./(x(:,3:end,:)-x(:,1:end-2,:));
dvdx(:,2:end-1,:) = (v(:,3:end,:)-v(:,1:end-2,:))./(x(:,3:end,:)-x(:,1:end-2,:));
dwdx(:,2:end-1,:) = (w(:,3:end,:)-w(:,1:end-2,:))./(x(:,3:end,:)-x(:,1:end-2,:));

dudy(2:end-1,:,:) = (u(3:end,:,:)-u(1:end-2,:,:))./(y(3:end,:,:)-y(1:end-2,:,:));
dvdy(2:end-1,:,:) = (v(3:end,:,:)-v(1:end-2,:,:))./(y(3:end,:,:)-y(1:end-2,:,:));
dwdy(2:end-1,:,:) = (w(3:end,:,:)-w(1:end-2,:,:))./(y(3:end,:,:)-y(1:end-2,:,:));

dudz(:,:,2:end-1) = (u(:,:,3:end)-u(:,:,1:end-2))./(z(:,:,3:end)-z(:,:,1:end-2));
dvdz(:,:,2:end-1) = (v(:,:,3:end)-v(:,:,1:end-2))./(z(:,:,3:end)-z(:,:,1:end-2));
dwdz(:,:,2:end-1) = (w(:,:,3:end)-w(:,:,1:end-2))./(z(:,:,3:end)-z(:,:,1:end-2));
    
%     %Loop through all coordinates in flow-field calculating the criteria
%     %for each point.
% for n = 2:S(1)-1
%     for m = 2:S(2)-1
%         for k = 2:S(3)-1
%                 %Calculates 3-by-3 Velocity Gradient Tensor
%             if flag2D   %Makes sure that incompressibility condition is satisfied for 2-D data.
%                 A = [dudx(n,m,k) dudy(n,m,k) dudz(n,m,k); dvdx(n,m,k) dvdy(n,m,k) dvdz(n,m,k); dwdx(n,m,k) dwdy(n,m,k) -dudx(n,m,k)-dvdy(n,m,k)];
%             else
%                 A = [dudx(n,m,k) dudy(n,m,k) dudz(n,m,k); dvdx(n,m,k) dvdy(n,m,k) dvdz(n,m,k); dwdx(n,m,k) dwdy(n,m,k) dwdz(n,m,k)];
%             end
% 
%             ls = eig(A); %eigenvalues of tensor
%             L(n,m,k) = max(unique(abs(imag(ls))));  %Swirling Strength
% 
%             Q(n,m,k) = -sum(sum(A.*A'))/2;    %Q-criterion
%             R = det(A);
%             D(n,m,k) = (Q(n,m,k)/3)^3 +(R/2)^2;  %Delta-criterion
% 
%             SR = (A +A')/2; %Strain rate tensor
%             OR = (A -A')/2; %Vorticity tensor
%             ls = sort(eig(SR^2 +OR^2));
%             L2(n,m,k) = ls(2);  %\lambda_2-criterion
%         end
%     end
% end

K = numel(u);
parfor k = 1:K
		%Calculates 3-by-3 Velocity Gradient Tensor
	if flag2D   %Makes sure that incompressibility condition is satisfied for 2-D data.
		A = [dudx(k) dudy(k) dudz(k); dvdx(k) dvdy(k) dvdz(k); dwdx(k) dwdy(k) -dudx(k)-dvdy(k)];
	else
		A = [dudx(k) dudy(k) dudz(k); dvdx(k) dvdy(k) dvdz(k); dwdx(k) dwdy(k) dwdz(k)];
	end

	ls = eig(A); %eigenvalues of tensor
	L(k) = max(unique(abs(imag(ls))));  %Swirling Strength
end

L = reshape(L,S);

%% OUTPUT Formatting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(dli==dlo)~=3 %Reorders the matrices to conform to input format
    cx = zeros(3,1);
    for n = 1:3
        cx(n) = find(dlo==dli(n));
    end
    
    L = permute(L,cx);
end 

if flag2D   %Removes erroneous data points in the instance of 2-D data
    L = L(:,:,2);
end

return

