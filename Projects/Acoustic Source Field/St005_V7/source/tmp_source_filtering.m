

L = size(phavg,3);

for n = 1:L
    tmp = phavg(:,:,n);
    
    tmp(1:37,1:2) = 0;
    tmp(1:70,end-1:end) = 0;
    tmp(1,end-4:end) = 0;
    tmp = moving_average2(tmp,1);
    pcolor(z/D,r/D,tmp/c/c);shading interp; colormap jet;colorbar;caxis([-2e7 2e7]);
    
    f = num2str(n);
    saveas(gcf,f,'fig');
    export_fig(f,'-png','-transparent','-r150');       
    
end