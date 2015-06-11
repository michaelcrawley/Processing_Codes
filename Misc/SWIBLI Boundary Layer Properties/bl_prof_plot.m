function bl_prof_plot(y,u,v_infty,M_infty,delta,theta)
%This function plots the boundary layer profile as normalized by Maise and
%McDonald (1968)

global To;
global r;
Tinf=To*(1+0.2*M_infty^2)^-1;
Ta=Tinf*(1+r*0.2*M_infty^2);
nu_w=-1.1555e-14*Ta^3+9.5728e-11*Ta^2+3.7604e-08*Ta-3.4484e-06;


%%%%%%%% Normalized Velocity Profile %%%%%%%%
eqn=strcat('C*(x*',num2str(delta/1000),'/',num2str(theta),')^(1/7)');
global fit
fit=ezfit(y./delta.*1000,u./v_infty,eqn);

figure
set(gcf,'Color','White')
plot(u./v_infty,y./delta.*1000,'+k');
title('Normalized Velocity Profile');
xlabel('u/u_\infty');
ylabel('y/\delta');
global dir_name;
set(gca,'XLim',[0,1.1]);
fn=strcat(dir_name,'\','normalized_vprofile');
saveFigure(gcf,fn,400);


%%%%%%%% u+/y+ Profile %%%%%%%%