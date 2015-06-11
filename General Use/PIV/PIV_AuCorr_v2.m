function [Cx,VCm] = PIV_AuCorr_v2(x,y,U,V,Um,Vm,L1,SW)

%L1 - File prefix
%SW - Potential core only logical

%Enter Um and Vm as zero to process single image

S = size(U);
qq = diff(sign(y(1,:)));
yi = (1:length(y(1,:))-1);
zr = yi(logical(qq));
CR = zr;

zUm = false;
if max(Um(:))==0
    zUm = true;
    Um = U;
    Vm = V;
    S(3) = 1;
end
if SW
    Q = Um(zr,:);  %grab centerline profile
    CI = max(find( Q > max(Q)*0.9));    %find index for end of potential core
    CI = round(CI*0.85);    %reduces core length to ensure we are in core region
    clear Q
    
    [C,I] = max(Um(CI,:));  %Finds peak velocity at end of potential core
    CL = min(find(Um(CI,:)>C/2));
    CL = I - 2*(I-CL);      %Sets top of potential core as center - 2*HWHM
%     CR = 2*I -CL;           %Sets bottom of potential core as center + 2*HWHM
else
    CI = length(x(:,1));
    CL = 1;
%     CR = S(2);
end

U = U(1:CI,CL:CR,:);
Um = Um(1:CI,CL:CR);
V = V(1:CI,CL:CR,:);
Vm = Vm(1:CI,CL:CR);

VC = []; UC = [];
if zUm
    Um = 0;
    Vm = 0;
end
for n = 1:S(3)
    progress(n,1,S(3),10);
    dU = U(:,:,n)'-Um';
    UC(:,:,n) = xcorr2(dU);
    UC(:,:,n) = UC(:,:,n)/max(max(UC(:,:,n)));
    
    dV = V(:,:,n)'-Vm';
    VC(:,:,n) = xcorr2(dV);
    VC(:,:,n) = VC(:,:,n)/max(max(VC(:,:,n)));
end
UCm = mean(UC,3);
UCs = std(UC,0,3);

VCm = mean(VC,3);
VCs = std(VC,0,3);

Cx = linspace(-x(CI,1),x(CI,1),2*CI-1);
Cy = linspace(-y(1,CL),y(1,CL),2*length(y(1,CR:CL))-1);

% L1 = 'm1';
save([cd filesep L1 '_Corr.mat'],'UC','VC','UCm','VCm','UCs','VCs','Cx','Cy')

% contourf(Cx,Cy,VCm)
% xlabel('Streamwise Shift (x/D)')
% ylabel('Radial Shift (y/D)')
% saveas(gcf,[cd L1 '_Corr.fig'])
% saveas(gcf,[cd L1 '_Corr.png'])