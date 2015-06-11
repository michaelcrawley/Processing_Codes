function [cylgrid,cylvec] = Cart2CylVectors(cartgrid,cartvec)
%Transforms Cartesian coordinates and vectors into Cylindrical
%Note that the standard definition, per Rachelle's LES simulation mesh, is
%that the jet axis dimension is labeled as 'x', not 'z'

    %Transform grid 
    [cylgrid.theta,cylgrid.r,cylgrid.z] = cart2pol(cartgrid.z,cartgrid.y,cartgrid.x);
    
    %Transform vectors
    cylvec.r = cos(cylgrid.theta).*cartvec.w + sin(cylgrid.theta).*cartvec.v;
    cylvec.theta = -sin(cylgrid.theta).*cartvec.w + cos(cylgrid.theta).*cartvec.v;
    cylvec.z = cartvec.u;
end