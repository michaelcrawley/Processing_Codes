% plot_times = [1 600 800 1000 6000];
% dx = 1; dt = 0.1;
plot_times = [1 1200 1600 2000 12000];
dx = 0.5; dt = 0.05;
dir = uigetdir([pwd,'\']);
cd(dir);
for i = 1:length(plot_times)
%    h=openfig(['density_i',num2str(plot_times(i))]);
%    saveas(gcf,['density_i',num2str(plot_times(i))],'png');
%    close(h);
%   
%    h=openfig(['u_i',num2str(plot_times(i))]);
%    saveas(gcf,['u_i',num2str(plot_times(i))],'png');
%    close(h);
%    
%    h=openfig(['v_i',num2str(plot_times(i))]);
%    saveas(gcf,['v_i',num2str(plot_times(i))],'png');
%    close(h);
   
   h=openfig(['p_i',num2str(plot_times(i))]);
   title(['Pressure contour for t = ',num2str(dt*plot_times(i)),', dx = ', num2str(dx),', dt = ', num2str(dt)]);
   saveas(gcf,['p_i',num2str(plot_times(i))],'png');
   saveas(gcf,['p_i',num2str(plot_times(i))],'png');
   close(h);   
end
cd ..