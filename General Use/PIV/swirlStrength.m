function L = swirlStrength(x,y,u,v)
%Calculates the swirling strength in velocity field using a second order
%accurate [(delta x)^2/3] velocity gradient tensor.  The complex component
%of the eigenvalues of this tensor is the swirling strength with units
%(1/seconds).
%
%Inputs are:    x - x-coordinate matrix
%               y - y-coordinate matrix
%               u - u-velocity matrix (total velocity - not fluctuating component)
%               v - v-velocity matrix

if x(2,1)-x(1,1)==0
    x = x'; y = y';
    u = u'; v = v';
end

S = size(u);
L = zeros(S(1),S(2));
    %calculates 2-by-2 tensor for each grid point
for n = 2:S(1)-1
    for m = 2:S(2)-1
            %Velocity Gradient Tensor
%         A = [(u(n,m)-u(n-1,m))/(x(n,m)-x(n-1,m)), (v(n,m)-v(n-1,m))/(x(n,m)-x(n-1,m));...
%             (u(n,m)-u(n,m-1))/(y(n,m)-y(n,m-1)), (v(n,m)-v(n,m-1))/(y(n,m)-y(n,m-1))];  %First order accurate    
        A = [(u(n+1,m)-u(n-1,m))/(x(n,m)-x(n-1,m))/2, (v(n+1,m)-v(n-1,m))/(x(n,m)-x(n-1,m))/2;... 
            (u(n,m+1)-u(n,m-1))/(y(n,m)-y(n,m-1))/2, (v(n,m+1)-v(n,m-1))/(y(n,m)-y(n,m-1))/2];  %Second order accurate
        
        Q = eig(A); %eigenvalues of tensor
        L(n,m) = abs(imag(Q(1)));
    end
end
