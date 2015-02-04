function [AO,AE,AW,AN,AS,S] = linkY(U,V,P,UR,VR,VL,VT,VB,X,Y,rho,mu)

%this function calculates the link coefficients for the Y-momentum equation
%given a guess for the velocity and pressure fields, and the mesh. Note
%that this is the specifically for a 2D staggered grid with the U-cell
%centered on the western boundary of the pressure-cell and the V-cell
%centered on the southern boundary. The code returns the link coefficients
%and source term for the Y momentum equation.
%
%Please note that this code calculates the coefficients as if for the X
%momentum equation, the switchs the directions at the end of the code. This
%was to ease the transfer from the x to the y momentum equation


%calling the X link coefficient generator
[AO,AN,AS,AE,AW,S] = linkX(V',U',P',VT',VB',VR',VL',UR',Y',X',rho',mu');


%transposing the matrices
AO = AO';
AE = AE';
AW = AW';
AN = AN';
AS = AS';
S = S';