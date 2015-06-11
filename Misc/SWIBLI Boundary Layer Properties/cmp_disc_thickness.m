function [cdisct,dr] = cmp_disc_thickness(u,y,Uinfty,Minfty)
%When given a vector of the velocities and vertical locations of the
%velocities (with y=0 being the surface) and the freestream velocity and
%Mach number, this function calculates the compressible displacement
%thickness of the boundary layer for which the data is provided. The
%formulas used to account for density variation can be found in Maise G.
%and McDonald, H. "Mixing Length and Kinematic Eddy Viscosity in a
%Compressible Boundary Layer." AIAA Journal Vol. 6:1, 1968. Please note
%that the value of r used in this paper is 1, while "Boundary Layer Theory"
%by Schlichting cites the value to be 0.89.

global r
dr(1)=(1+r*0.2*Minfty.^2*(1-u(1).^2/Uinfty.^2)).^-1;
cdisct=0;
for i=2:1:length(u)
    %dr is the density ratio
    dr(i)=(1+r*0.2*Minfty.^2*(1-u(i).^2/Uinfty.^2)).^-1;
    cdisct=cdisct+(1-dr(i)*u(i)/Uinfty)*(y(i)-y(i-1));
end