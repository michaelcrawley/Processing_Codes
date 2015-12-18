function delta=bl_prop(u_infty,M_infty,x,percent,u,y)
%When given a vector of normalized velocities, and the vertical coordinates of
%those velocities and the freestream velocity this function will calculate
%the boundary layer thickness, and
%the n for the best-fit log-law velocity profile for the boundary layer, it
%will plot the original points, the best-fit log-law profile, and a 1/7
%log-law profile in normalized form. It will then calculate boundary layer
%properties such as the displacement thickness, momentum thickness, shape
%factor, and Reynolds number based on momentum thickness.
%
%This function takes 4 arguments u_infty-the freestream velocity;
%percent-the percent of the freestream velocity that is used as the cutoff
%for the boundary layer; u-a vector of velocities (in m/s); y-a vector of
%vertical locations of the velocities
%
%x is the streamwise location of the boundary layer measurement relative
%tot he end of the nozzle (that is the point at which it begins to develop
%under a constant pressure

global logf

u_n=u./u_infty;

i=1;
while(percent/100>u_n(i) && i<length(u_n))
    i=i+1;
end
delta=y(i-1);

y_n=y./delta; %mm

%calculation incompressible shape factor
idist=inc_disc_thickness(u,y./1000,u_infty);
imomt=inc_mom_thickness(u,y./1000,u_infty);
Hinc=idist/imomt;
fprintf('The incompressible displacement thickness is %.3f mm\n',idist*1000);
fprintf('The incompressible momentum thickness is %.3f mm\n',imomt*1000);
fprintf('The incompressible shape factor is %.3f\n\n',Hinc);

fprintf(logf,'The incompressible displacement thickness is %.3f mm\n',idist*1000);
fprintf(logf,'The incompressible momentum thickness is %.3f mm\n',imomt*1000);
fprintf(logf,'The incompressible shape factor is %.3f\n\n',Hinc);

%calculation compressible shape factor
[cdist,dr]=cmp_disc_thickness(u,y./1000,u_infty,M_infty);
[cmomt,dr]=cmp_mom_thickness(u,y./1000,u_infty,M_infty);
Hcmp=cdist/cmomt;
fprintf('The compressible displacement thickness is %.3f mm\n',cdist*1000);
fprintf('The compressible momentum thickness is %.3f mm\n',cmomt*1000);
fprintf('The compressible shape factor is %.3f\n\n',Hcmp);

fprintf(logf,'The compressible displacement thickness is %.3f mm\n',cdist*1000);
fprintf(logf,'The compressible momentum thickness is %.3f mm\n',cmomt*1000);
fprintf(logf,'The compressible shape factor is %.3f\n\n',Hcmp);


bl_prof_plot(y./1000,u,u_infty,M_infty,delta,cmomt);


%calculation reynolds number based on the momentum thickness
global To
global Po
T=To/(1+0.2*M_infty^2);
%viscosity is calculated using southerland's law
Tos=524.07; %degR
muos=0.01827; %centipoise
C=120; %southerland's constant
TR=((9/5*(T+273.15))+32)-459.67;
a=0.555*Tos+C;
b=0.555*TR+C;
mu=0.001*muos*(a/b)*(TR/Tos)^(3/2);

%calculating the freestream flow density
rhoo=Po/287/To; %kg/m^3
rho=rhoo*(1+0.2*M_infty^2)^-2.5;

Re_theta=rho*u_infty*cmomt/mu;
Re_x=rho*u_infty*x/mu/1000;
fprintf('The boundary layer thickness is %.3f mm\n',delta);
fprintf('The Reynolds number based on the momentum thickness is %.f\n',Re_theta);
fprintf('The Reynolds number based on the streamwise location is %.3e\n\n',Re_x);

fprintf(logf,'The boundary layer thickness is %.3f mm\n',delta);
fprintf(logf,'The Reynolds number based on the momentum thickness is %.f\n',Re_theta);
fprintf(logf,'The Reynolds number based on the streamwise location is %.3e\n\n',Re_x);

%plotting the density ratio profile
figure
set(gcf,'Color','White');
plot(dr,y_n,'-k');
title('Density Ratio Profile');
xlabel('\rho/\rho_\infty');
ylabel('y/\delta');
set(gca,'YLim',[0,1]);
global dir_name
fn=strcat(dir_name,'\','density_ratio_profile');
saveFigure(gcf,fn,400);