
%- Load PIV processed data set
%- Run PIV_CorrMap to get conditionally averaged images
%- Run VortexID to get structure information
%- Find peaks in swirl strength to measure structure: size, strength, and
%convective velocity

D = 0.0254; %m - jet diameter

%%%%% For low frequency forcing
% [CondAvg.Matches,CondAvg.CorrCoeff,CondAvg.StructSpace] = PIV_CorrMap_v3(x,y,Um,V,Vm,1); clear lastfit;

%%%%% High frequency forcing
xl = [1 find(x(:,1) > 3,1,'first')];
yl = [find(y(1,:) > 0.6,1,'last') find(y(1,:) < -0.6,1,'first')];
[CondAvg.Matches,CondAvg.CorrCoeff,CondAvg.StructSpace] = PIV_CorrMap_v3(x,y,Um,V,Vm,[xl yl]); clear lastfit;

L = zeros([size(Um) size(CondAvg.Matches,1)]); CondAvg.L = L; CondAvg.Q = L; 
CondAvg.L2 = L; CondAvg.U = L; CondAvg.V = L; clear L;
for n = 1:size(CondAvg.Matches,1)
    M = CondAvg.Matches(n,:);
    CondAvg.V(:,:,n) = mean(V(:,:,M(M>0)),3);   %creates conditionally averaged whole images
    CondAvg.U(:,:,n) = mean(U(:,:,M(M>0)),3);
    CondAvg.TKE(:,:,n) = (std(U(:,:,M(M>0)),0,3).^2 +std(V(:,:,M(M>0)),0,3).^2)/Uj^2;
    [Dl,CondAvg.L(:,:,n),CondAvg.Q(:,:,n),CondAvg.L2(:,:,n)] = VortexID(x*D,y*D,CondAvg.U(:,:,n),CondAvg.V(:,:,n)); 
    clear M Dl
end
clear U V

xl = [1 find(x(:,1) > 8,1,'first')];    %x-limits of structure region
yl = [find(y(1,:) > 2,1,'last') find(y(1,:) < -2,1,'first')];   %y-limits of structure region
Uc = CondAvg.U(xl(1):xl(2),yl(1):yl(2),:);  %Temp variable holding u-velocity in structure region
Q = CondAvg.L(xl(1):xl(2),yl(1):yl(2),:);   %Temp variable holding swirling strength in structure region
cl = [mean(Q(:))+std(Q(:))*2 max(Q(:))];  %color limits

done = false;   
while ~done
    for n = 1:size(CondAvg.Matches,1)
        pcolor(Q(:,:,n)'); shading flat; caxis([0 cl(2)]); colorbar;
        hold on
        [CM,h] = contour(Q(:,:,n)',[-1 cl(1)]); set(h,'LineColor','w','LineWidth',2); pause(1);
        hold off
    end
    done = input('Is it done (1/0): ');
    if ~done
        cl(1) = input(['What threshold (last = ' num2str(cl(1)) '): ']);
        figure(gcf);
    end
end
close
clear n done CM h
I = Q >= cl(1);  %Indices of vortices
tp1 = zeros(size(Uc(:,:,1))); tp2 = tp1; TP1 = tp1;
for n = 1:size(CondAvg.Matches,1)
    tmp = Uc(:,:,n);
    tp1(I(:,:,n)) = tp1(I(:,:,n))+tmp(I(:,:,n));
    tp2(I(:,:,n)) = tp2(I(:,:,n))+1;
    
    tmp = Q(:,:,n);
    TP1(I(:,:,n)) = TP1(I(:,:,n))+tmp(I(:,:,n));
end
tp3 = tp1./tp2; TP3 = TP1./tp2; clear tp1 tp2 tmp TP1;
tp4 = zeros(xl(2),1); TP4 = tp4;
for n = 1:xl(2)
    tp4(n) = mean(tp3(n,~isnan(tp3(n,:))));
    TP4(n) = mean(TP3(n,~isnan(TP3(n,:))));
end

CondAvg.AvgL = TP4;  %Average swirling Strength (1/s)
CondAvg.AvgUc = tp4;    %Average convective velocity (m/s)
clear n I Q Uc tp3 tp4 TP3 TP4


% save('StructInfo.mat','x','y','CorrCoeff','C','CUavg','CVavg','L','Q','L2','SS','Uc','Vn','Vm','Vp','Vu','Vx','Vy');



%%%% PLOTTING %%%%%%%%%%%%%
% % % contourf(x,y,L,12)
% % % axis([min(x(:)) 8 -2 2])
% % % xlabel('x/D')
% % % ylabel('y/D')
% % % grid on
% % % colormap cool
% % % colorbar
% % % title({'Flow Structures: T_o/T_a \approx 1.5, m = \pm1',...
% % %     'Uc = 187 m/s, Structure Spacing = 1.62 x/D',...
% % %     'Vortex Peak Swirling Strength = 270 1/sec'})
% % % saveas(gcf,'Vm11_Swirl.fig')
% % % saveas(gcf,'Vm11_Swirl.png')
% % % 
% % %     %For bad averaging
% % % save('Base_StructInfo.mat','x','y','C','CUavg','CVavg','L');
% % % contourf(x,y,L,12)
% % % axis([min(x(:)) 8 -2 2])
% % % xlabel('x/D')
% % % ylabel('y/D')
% % % grid on
% % % colormap cool
% % % colorbar
% % % title({'Flow Structures: T_o/T_a = 1.5, Baseline',...
% % %     'No Acceptable Result From Conditional Averaging',...
% % %     'Vortex Peak Swirling Strength \approx 146 1/sec'})
% % % saveas(gcf,'Base_Swirl.fig')
% % % saveas(gcf,'Base_Swirl.png')
% % % 
% % % 
% % % %%% For Figures with Galilean streamlines
% % % contourf(x,y,L,12)
% % % axis([min(x(:)) 9 -2 2])
% % % xlabel('x/D')
% % % ylabel('y/D')
% % % grid on
% % % colormap summer
% % % colorbar
% % % 
% % % Res = 0.1;
% % % xm = find(x(:,1) > 9,1,'first')-1;
% % % hold on
% % % for qy = -1.2:Res:-.2
% % %     h = streamline(x(1:xm,:)',y(1:xm,:)',CUavg(1:xm,:)'-Uc,CVavg(1:xm,:)',8.5,qy);
% % %     set(h,'Color','black')
% % % end
% % % for qy = -0.5:Res:0.5
% % %     h = streamline(x(1:xm,:)',y(1:xm,:)',CUavg(1:xm,:)'-Uc,CVavg(1:xm,:)',x(1,1),qy);
% % %     set(h,'Color','black')
% % % end
% % % for qy = 0.5:Res:1.2
% % %     h = streamline(x(1:xm,:)',y(1:xm,:)',CUavg(1:xm,:)'-Uc,CVavg(1:xm,:)',8.5,qy);
% % %     set(h,'Color','black')
% % % end
% % % hold off
% % % clickStreamline(x(1:xm,:),y(1:xm,:),CUavg(1:xm,:)-Uc,CVavg(1:xm,:),8.5,'black')
% % % axis([min(x(:)) 8 -1.5 1.5])
% % % 
% % % saveas(gcf,'Base_Swirl.fig')
