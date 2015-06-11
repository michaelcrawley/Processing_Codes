%Throwaway function for creating animated gif


for i = a:3:b
    h = figure;surf(100*x, 100*x',Tn(:,:,i)-Ta(:,:,i));shading interp;title(['time: ' num2str(t(i))]);
    zlim([min(min(Tn(:,:,a)-Ta(:,:,a))) max(max(Tn(:,:,a)-Ta(:,:,a)))]);
    gif_add_frame(h,filename,24);
    close(h);
end