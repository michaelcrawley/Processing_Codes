%BL height vs. z-position
figure(8)
set(gcf,'color','w')
%define z plane limits
for l = 1 : size(z{i},1)
    if z{i}(l,1)<-lm & z{i}(l+1,1) > -lm
        mnz=l;
        break
    end
end
subplot(4,2,i)
plot(z{i}(mnz:mnz+np1(i)-1,1),ybl1{i});
hold all
% plot(z1{1}(1:np12,1),ybl12);
plot(z{i}(mnz:mnz+np1(i)-1,1),tot_disct{i}(mnz:mnz+np1(i)-1));
plot(z{i}(mnz:mnz+np1(i)-1,1),tot_momt{i}(mnz:mnz+np1(i)-1));
plot(z{i}(mnz:mnz+np1(i)-1,1),tot_H{i}(mnz:mnz+np1(i)-1));
locp=[ 22 47 52 57 62 67 72 77];
title(['Boundary-Layer properties vs. Z-Position, x = ', num2str(locp(i)), ' mm'],'fontname','times new roman','fontsize',12)
ylim([0,12]);
ylabel('\delta, \delta*, \theta, H','fontname','times new roman','fontsize',12)
xlabel('z (mm)','fontname','times new roman','fontsize',12)
leg_e{1}='BL height (Max RMS)';
%leg_e{2}='BL1 height (Mean)';
leg_e{2}='BL displacement thickness';
leg_e{3}='BL momentum thickness';
leg_e{4}='BL shape factor';
legend(leg_e);
% file=(['fig\Bl_properties_x',num2str(locp(i))]);
% saveas(gcf,file,'fig')
% file=(['bmp\Bl_properties_x',num2str(locp(i))]);
% saveas(gcf,file,'bmp')
% clf