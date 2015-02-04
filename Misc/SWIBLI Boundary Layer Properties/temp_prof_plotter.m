%Temporary plotting file
close all
clear
clc

M_infty=2.045;
v_infty=514.3;
delta=5.099;
theta=0.62;
r=0.896;

load BL_prof.mat
y=vert./1000;
u=vel;

% u=v_infty.*[0,.4078,.6531,.8449,.9264,.9607,.981,.9901];
% y=[0,.4812,1.283,2.085,2.887,3.689,4.491,5.293]./1000; %mm


% figure
% set(gcf,'Color','White');
% set(gcf,'Position',[149,436,1126,560]);
% subplot(1,2,1)
% plot(y./delta.*1000,u./v_infty,'+')
% set(gca,'Position',[0.075,0.11,0.38,0.815]);
% set(gca,'XLim',[0,1.2])
% subplot(1,2,2)
% semilogx(y./delta.*1000,u./v_infty,'+')
% set(gca,'Position',[0.54,0.11,0.38,0.815]);
% set(gca,'XLim',[0,1.2])
% removeGraySpace(gcf)



%%%%%%%%%%%%%%%%%%%   U*, Maise and McDonald   %%%%%%%%%%%%%%%%%%%%%%%%%


%define these (assume adiabatic)
To=298.15;%252.3; %K
Tinf=To*(1+0.2*M_infty^2)^-1;
Tw=Tinf*(1+r*0.2*M_infty^2);
nu_w=-1.1555e-14*Tw^3+9.5728e-11*Tw^2+3.7604e-08*Tw-3.4484e-06;
%calc dudy
maxi=1;
dudy=0;
i=9;
%for i=1:maxi
    dudy=(-u(i+2)+4*u(i+1)-3*u(i))/(y(i+2)-y(i));
%end
%dudy=dudy./maxi;

u_tau=14.33;%44.37;%sqrt(nu_w*dudy);
% A=sqrt(0.2*M_infty^2/(Tw/Tinf));
% B=(1+0.2*M_infty^2)/(Tw/Tinf)-1;
% u_inf_star=v_infty/A*asin((2*A^2-B)/sqrt(B^2+4*A^2));
% u_star=v_infty/A.*asin((2*A^2.*(u./v_infty)-B)./sqrt(B^2+4*A^2));
% u_normalized=(u_inf_star-u_star)./u_tau;
% 
% y_normalized=y./delta.*1000;
% 
% %model
% y_mod_n=[0.01:0.01:1];
% u_mod_n=-2.5.*log(y_mod_n)+1.25.*(1+cos(3.1415926.*y_mod_n));
% 
% 
% % %portion over which to fit
% % s=1;
% % e=length(u);
% % 
% % 
% % a=1;
% % b=100;
% % step=1;
% % while step>=0.001
% %     sum_min=inf;
% %     min=1;
% %     for ut=a:step:b
% %         sum=0;
% %         for i=s:e
% %             sum=sum+((u_inf_star-u_star(i))/ut-(-2.5*log(y_normalized(i))+1.25*(1+cos(3.1415926*y_normalized(i)))))^2;
% %         end
% %         sum=sqrt(sum);
% % 
% %         if sum<sum_min
% %             sum_min=sum;
% %             min=ut;
% %         end
% %     end
% %     a=min-step;
% %     b=min+step;
% %     step=step/10;
% % end
% % 
% % u_tau=min;
% u_normalized=(u_inf_star-u_star)./u_tau;
% fprintf('\n>> \n\nU_tau=%.3f\n\n',u_tau);
% 
% 
% 
% figure
% set(gcf,'Color','White');
% semilogx(y_normalized(1:55),u_normalized(1:55),'+k',y_mod_n,u_mod_n,'-k');
% title('Experimental and Generalized Literature Boundary Layer Profiles');
% xlabel('y/\delta');
% ylabel('(u^*_\infty-u^*)/u_\tau');
% legend('Experimental','Literature Profile');
% set(legend,'Location','NorthEast');
% set(gca,'XLim',[0.01,1]);
% set(gca,'XTick',[0.01,0.02,0.05,0.1,0.2,0.5,1.0]);
% % set(gca,'YLim',[0,16]);
% % set(gca,'YTick',[0:2:16]);




%%%%%%%%%%%%%%%%%%    U+, y+    %%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum_min=inf;
% for ut=1:1:100
%     sum=0;
%     for i=7:33%1:1:length(u)
%         sum=sum+(u(i)./ut-(5.2+0.41.^-1.*log(y(i)./(nu_w/u_tau))))^2;
%     end
%     sum=sqrt(sum);
%     if sum<sum_min
%         u_tau=ut;
%         sum_min=sum;
%     end
% end


fprintf('\n>> \n\nU_tau=%.2f m/s\n\n',u_tau);

u_n=u./u_tau;
y_n=y.*(u_tau/nu_w);

y_mod_n=y_n;
u_mod_n=5.2+0.41^-1.*log(y_mod_n);%13,0.36,26

% fit=ezfit(y_n(8:32),u_n(8:32),'C+1/k*log(x/r)');
% 
% yfit=[1000:10:4000];
% ufit=fit.m(1)+fit.m(2)^-1.*log(yfit./fit.m(3));

figure
set(gcf,'Color','White')
semilogx(y_n,u_n,'+k',y_mod_n,u_mod_n,'-k');%,yfit,ufit,'--k');
%set(gca,'XLim',[10,8000]);
title('U^+ vs. y^+:Experimental and Literature Profiles');
xlabel('y^+=y/\delta_\nu');
ylabel('u^+=u/u_\tau');
legend('Experimental','Literature Profile');%,'Best-Fit Profile');
set(legend,'Location','NorthWest');
set(gca,'XLim',[100,20000]);
set(gca,'YLim',[0,35]);
grid on






% ym=[0:0.01:14];
% um=0.683.*(ym).^(1/7);
% ym=[ym,[14:0.01:22]];
% um=[um,ones(1,length(14:0.01:22))];
% figure
% plot(ym,um,'--k',y./theta.*1000,u./v_infty,'+');