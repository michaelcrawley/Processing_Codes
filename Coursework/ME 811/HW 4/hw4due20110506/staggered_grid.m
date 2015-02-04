function [Xp,Yp,Xu,Yu,Xv,Yv] = staggered_grid(Lx,Ly,M,N)

%staggered_grid generates the u-velocity, v-velocity, and pressure grids
%for a 2D staggered mesh for use for the SIMPLE algorithm. It takes the
%physical x and y dimensions of the grid, and the number of volumes in the
%pressure mesh in the x and y directions:
%
%[Xu,Yu,Xv,Yv,Xp,Yp] = staggered_grid(Lx,Ly,M,N)
%
%Note that the outputs to this are matrices: the y-direction varies with
%the row index, and the x-direction with the column index.


%initializing the variables
xu = zeros(1,M);
yu = zeros(N,1);
xv = zeros(1,M);
yv = zeros(N,1);
xp = zeros(1,M);
yp = zeros(N,1);


%generating the pressure grid
dx = Lx/M;
dy = Ly/N;

xp = [dx/2:dx:Lx-dx/2];
yp = [dy/2:dy:Ly-dy/2];


%generating the velocity grids
xu = xp - dx/2;
yu = yp;

xv = xp;
yv = yp - dy/2;


%generating the grid matrices
[Xp,Yp] = meshgrid(xp,yp);
[Xu,Yu] = meshgrid(xu,yu);
[Xv,Yv] = meshgrid(xv,yv);