function [xc,AC] = PIV_AutoCorr(x,y,Um,V,Vm,SW)
%Calculates the V-fluctuating-component auto-correlation in the streamwise
%direction.
%
%Inputs:    x - x-coordinate matrix
%           y - y-coordinate matrix
%           V - v-velocity matrix - streamwise-by-radial-by-Image
%           Um - u-velocity average image
%           Vm - v-velocity average image
%           SW - switch parameter 
%               - 1==Correlate potential core only
%               - 0==Correlate whole image
%               - else == Correlate over this distance
%
%Outputs:   xc - displacement coordinate of correlation
%           AC - auto-correlation



S = size(V);
qq = diff(sign(y(1,:)));
yi = (1:length(y(1,:))-1);
zr = yi(logical(qq));   %Finds index of jet centerline
clear qq yi
switch length(SW)
    case 1
        CO = 1;
        if SW==1
            Q = Um(:,zr);  %grab centerline profile
            CI = find( Q > max(Q)*0.9,1,'last');    %find index for end of potential core
            CI = round(CI*0.85);    %reduces core length to ensure we are in core region
            disp(['Correlating in potential core region x/D < ' num2str(x(CI,1))])
            clear Q

            [C,I] = max(Um(CI,:));  %Finds peak velocity at end of potential core
            CL = find(Um(CI,:)>C/2,1,'first');
            CL = I - 2*(I-CL);      %Sets top of potential core as center - 2*HWHM
            CR = 2*I -CL;           %Sets bottom of potential core as center + 2*HWHM
        elseif SW==0
            CI = length(x(:,1));
            CL = 1;
            CR = S(2);
        else
            CI = find(x>SW,1,'first');
            [C,I] = max(Um(CI,:));  %Finds peak velocity at end of potential core
            CL = find(Um(CI,:)>C/2,1,'first');
            CL = I - 2*(I-CL);      %Sets top of potential core as center - 2*HWHM
            CR = 2*I -CL;           %Sets bottom of potential core as center + 2*HWHM
        end
    case 2
        CO = min(SW); CI = max(SW); CL = [0 0]; CR = CL;
        for n = 1:2
            CI = find(x>SW(n),1,'first');
            [C,I] = max(Um(CI,:));  %Finds peak velocity at end of potential core
            CL(n) = find(Um(CI,:)>C/2,1,'first');
            CL(n) = I - 2*(I-CL(n));      %Sets top of potential core as center - 2*HWHM
            CR(n) = 2*I -CL;           %Sets bottom of potential core as center + 2*HWHM
        end
        CL = min(CL);
        CR = max(CR);
    case 4
        CO = SW(1); CI = SW(2); CL = SW(3); CR = SW(4);
    otherwise
    error('bad switch')
end
    %crops velocity maps if SW~=0 so that correlation only looks at jet
    %potential core
V = V(CO:CI,CL:CR,:);
Vm = Vm(CO:CI,CL:CR);



AC = zeros(size(Vm)*2-1);
% for n = 1:S(3)-1
%     progress(n,1,S(3)-1,10);
%     for m = n+1:S(3)
%         c = c+1;
%         C = xcorr2(V(:,:,n)-Vm,V(:,:,m)-Vm);
%         XC = XC+C/max(C(:));
%     end
% end
% XC = XC/c;
for n = 1:S(3)
    C = xcorr2(V(:,:,n)-Vm,V(:,:,n)-Vm);
    AC = AC+C/max(C(:));
end
AC = AC/S(3);
AC = AC(:,size(Vm,2));

dx = abs(mean(diff(x(:,1))));
xc = (-dx*(size(Vm,1)-1):dx:dx*(size(Vm,1)-1))';




