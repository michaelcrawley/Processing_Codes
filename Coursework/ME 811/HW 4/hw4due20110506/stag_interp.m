function [Uint,Vint] = stag_interp(U,V,UR,VT)

%stag_interp interpolates the velocity values onto the pressure grid for a
%standard staggered mesh when given the velocity field, and the values of
%the right-hand and top faces of the domain for the u and v velocities
%respectively. The interpolated velocities are returned as two matrices

%parsing the inputs
M = size(U,2);
N = size(U,1);
if size(UR,1) < size(UR,2)
    UR = UR';
end
if size(VT,1) > size(VT,2)
    VT = VT';
end

%assembling the matrices for interpolation
U = [U,UR];
V = [V;VT];

%intializing the interpolated matrices and performing the interpolation
Uint = zeros(N,M);
Vint = zeros(N,M);
for i = 1:M
    for j = 1:N
        Uint(j,i) = (U(j,i) + U(j,i+1))/2;
        Vint(j,i) = (V(j,i) + V(j+1,i))/2;        
    end
end