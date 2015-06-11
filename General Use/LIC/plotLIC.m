clear variables

fd = 'F:\Research\Mach1.3\PIV\Processed\20080930';
sf = {'m0_St013','m0_St026','m0_St033','m0_St052','m0_St105',...
    'm1_St013','m1_St026','m1_St033','m1_St052','m1_St105',...
    'vm11_St013','vm11_St026','vm11_St033','vm11_St052','vm11_St105'};
sn = {'m = 0, St_{DF} = 0.13','m = 0, St_{DF} = 0.26','m = 0, St_{DF} = 0.33','m = 0, St_{DF} = 0.52','m = 0, St_{DF} = 1.05',...
    'm = 1, St_{DF} = 0.13','m = 1, St_{DF} = 0.26','m = 1, St_{DF} = 0.33','m = 1, St_{DF} = 0.52','m = 1, St_{DF} = 1.05',...
    'm = \pm1, St_{DF} = 0.13','m = \pm1, St_{DF} = 0.26','m = \pm1, St_{DF} = 0.33','m = \pm1, St_{DF} = 0.52','m = \pm1, St_{DF} = 1.05'};

load([fd '\Distilled.mat'],'Ucfit');

for n = 1:length(sf)
    load([fd '\Condensed.mat'],sf{n})
    eval(['x = ' sf{n} '.x;']); eval(['y = ' sf{n} '.y;']);
    yl(1) = find(y(1,:)<1.5,1,'first')-1;
    yl(2) = find(y(1,:)<-1.5,1,'first');
    xl(1) = find(x(:,1)>0.6,1,'first')-1;
    xl(2) = find(x(:,1)>8,1,'first');
    
    x = x(xl(1):xl(2),yl(1):yl(2));
    y = y(xl(1):xl(2),yl(1):yl(2));
    eval(['u = ' sf{n} '.CondAvg.U(xl(1):xl(2),yl(1):yl(2),:);']); N = size(u,3);
    eval(['vG = ' sf{n} '.CondAvg.V(xl(1):xl(2),yl(1):yl(2),:);']);
    uG = u - repmat(Ucfit(xl(1):xl(2)),[1 size(u,2) N]);
    eval(['clear ' sf{n}])

    ph = (0:360/N:(N-1)/N*360);
    for m = 1:N
        [X,Y,SL(:,:,m)] = myStreamlines(x,y,uG(:,:,m),vG(:,:,m));
%         pcolor(X,Y,double(SL(:,:,m))); shading flat; colormap gray;
%         axis([0.6 8 -1.5 1.5])
%         xlabel('x/D')
%         ylabel('y/D')
        CL(2) = max(u(:)); CL(1) = CL(2)/4;
        NC =10; CM = makeLine([0 0 1],[0 1 0],ceil(NC/2)); CM = [CM(1:end-1,:); makeLine([0 1 0],[1 0 0],ceil(NC/2))];

        QS = SL(:,:,m);
        L = interp2(u(:,:,m)/CL(2),4);
        L = round(L*size(CM,1));
        I = L==0;
        Lrgb = ind2rgb(L,CM); clear L CM NC;
        Lycbcr = rgb2ycbcr(Lrgb); clear Lrgb;
        Lycbcr(:,:,1) = QS;
        Lrgb2 = ycbcr2rgb(Lycbcr); clear Lycbcr;
        tmp = Lrgb2(:,:,1); tmp(I) = QS(I); Lrgb2(:,:,1) = tmp;
        tmp = Lrgb2(:,:,2); tmp(I) = QS(I); Lrgb2(:,:,2) = tmp;
        tmp = Lrgb2(:,:,3); tmp(I) = QS(I); Lrgb2(:,:,3) = tmp; clear tmp QS;

        P = permute(Lrgb2,[2 1 3]); clear Lrgb2;
        image(X(:,1),Y(1,:),P)
        set(gca,'YDir','normal')
        axis([0.6 8 -1.5 1.5])
        xlabel('x/D')
        ylabel('y/D')
        title(['GS: T_o/T_a = 1.5, ' sn{n} ', \phi = ' num2str(ph(m))])
        saveFigure_v2(gcf,[fd '\GS_' sf{n} '_phi' num2str(ph(m))],600)
        close;
    end
    
    save([fd '\GS_' sf{n} '.mat'],'X','Y','SL','x','y','uG','vG','Ucfit');
end





% axis([1 8 -1.5 1.5])
% colormap(CM);
% colorbar;
% CT = round((CL(2)-CL(1))/NC*((1:NC)-0.5));  YTL = cell(1,NC);
% for n = 1:NC
%     YTL{n} = num2str(CT(n));
% end
% colorbar('YTick',(1:NC),'YTickLabel',YTL); clear P CM CT YTL CL;

% xlabel('x/D')
% ylabel('y/D')
