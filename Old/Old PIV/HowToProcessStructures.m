
%- Load PIV processed data set
%- Run PIV_CorrMap to get conditionally averaged images
%- Run swirlStrength to get structure information
%- Find peaks in swirl strength to measure structure: size, strength, and
%convective velocity

% clear TKE Tj To Uj CL
% [CUavg,CVavg,C,CorrCoeff] = PIV_CorrMap(x,y,U,Um,V,Vm,1);
[CondAvg.Matches,CondAvg.CorrCoeff,CondAvg.StructSpace] = PIV_CorrMap_v3(x,y,Um,V,Vm,1);
L = zeros([size(Um) size(CondAvg.Matches,1)]); CondAvg.L = L; CondAvg.Q = L; 
CondAvg.L2 = L; CondAvg.U = L; CondAvg.V = L; clear L;
for n = 1:size(CondAvg.Matches,1)
    M = CondAvg.Matches(n,:);
    CondAvg.V(:,:,n) = mean(V(:,:,M(M>0)),3);   %creates conditionally averaged whole images
    CondAvg.U(:,:,n) = mean(U(:,:,M(M>0)),3);
    [D,CondAvg.L(:,:,n),CondAvg.Q(:,:,n),CondAvg.L2(:,:,n)] = VortexID(x,y,CondAvg.U(:,:,n),CondAvg.V(:,:,n)); 
    clear M D
end
% imagesc(L')
% colormap jet
%Fill in Vn and Vm with vortex peak locations
%Vn = [];   %streamwise axis
%Vm = [];   %radial axis

Vn = Vn(:);
Vm = Vm(:);
Vn2 = [Vn Vn+1 Vn-1 Vn Vn Vn+1 Vn-1 Vn+1 Vn-1]; %averages 3 by 3 centered on peak
Vm2 = [Vm Vm Vm Vm+1 Vm-1 Vm+1 Vm+1 Vm-1 Vm-1];
Vp = mean(L(sub2ind(size(L),Vn2,Vm2)),2); %Vortex peak swirling strength (1/s)
Vu = mean(CUavg(sub2ind(size(L),Vn2,Vm2)),2); %Vortex peak u-velocity (m/s)
Vx = x(sub2ind(size(L),Vn,Vm));     %Vortex peak x-coordinate (x/D)
Vy = y(sub2ind(size(L),Vn,Vm));     %Vortex peak y-coordinate (y/D)

Uc = mean(Vu);    %Convective velocity (m/s)
mVp = mean(Vp);     %mean peak swirling strength (1/s)
SS = mean([diff(Vx(Vy > 0)); diff(Vx(Vy < 0))]); %Vortex streamwise spacing (x/D)

save('StructInfo.mat','x','y','CorrCoeff','C','CUavg','CVavg','L','Q','L2','SS','Uc','Vn','Vm','Vp','Vu','Vx','Vy');



%%%% PLOTTING %%%%%%%%%%%%%
contourf(x,y,L,12)
axis([min(x(:)) 8 -2 2])
xlabel('x/D')
ylabel('y/D')
grid on
colormap cool
colorbar
title({'Flow Structures: T_o/T_a \approx 1.5, m = \pm1',...
    'Uc = 187 m/s, Structure Spacing = 1.62 x/D',...
    'Vortex Peak Swirling Strength = 270 1/sec'})
saveas(gcf,'Vm11_Swirl.fig')
saveas(gcf,'Vm11_Swirl.png')

    %For bad averaging
save('Base_StructInfo.mat','x','y','C','CUavg','CVavg','L');
contourf(x,y,L,12)
axis([min(x(:)) 8 -2 2])
xlabel('x/D')
ylabel('y/D')
grid on
colormap cool
colorbar
title({'Flow Structures: T_o/T_a = 1.5, Baseline',...
    'No Acceptable Result From Conditional Averaging',...
    'Vortex Peak Swirling Strength \approx 146 1/sec'})
saveas(gcf,'Base_Swirl.fig')
saveas(gcf,'Base_Swirl.png')


%%% For Figures with Galilean streamlines
contourf(x,y,L,12)
axis([min(x(:)) 9 -2 2])
xlabel('x/D')
ylabel('y/D')
grid on
colormap summer
colorbar

Res = 0.1;
xm = find(x(:,1) > 9,1,'first')-1;
hold on
for qy = -1.2:Res:-.2
    h = streamline(x(1:xm,:)',y(1:xm,:)',CUavg(1:xm,:)'-Uc,CVavg(1:xm,:)',8.5,qy);
    set(h,'Color','black')
end
for qy = -0.5:Res:0.5
    h = streamline(x(1:xm,:)',y(1:xm,:)',CUavg(1:xm,:)'-Uc,CVavg(1:xm,:)',x(1,1),qy);
    set(h,'Color','black')
end
for qy = 0.5:Res:1.2
    h = streamline(x(1:xm,:)',y(1:xm,:)',CUavg(1:xm,:)'-Uc,CVavg(1:xm,:)',8.5,qy);
    set(h,'Color','black')
end
hold off
clickStreamline(x(1:xm,:),y(1:xm,:),CUavg(1:xm,:)-Uc,CVavg(1:xm,:),8.5,'black')
axis([min(x(:)) 8 -1.5 1.5])

saveas(gcf,'Base_Swirl.fig')
