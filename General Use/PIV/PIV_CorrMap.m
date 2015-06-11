function [CUavg,CVavg,C,CorrCoeff] = PIV_CorrMap(x,y,U,Um,V,Vm,SW)
%Calculates the V-fluctuating-component correlation map to create a
%conditionally averaged set of images which are 10% of the total set.
%
%Inputs:    x - x-coordinate matrix
%           y - y-coordinate matrix
%           U - u-velocity matrix - streamwise-by-radial-by-Image
%           V - v-velocity matrix - streamwise-by-radial-by-Image
%           Um - u-velocity average image
%           Vm - v-velocity average image
%           SW - switch parameter 
%               - 1==Correlate potential core only
%               - 0==Correlate whole image
%
%Outputs:   CUavg - Conditionally averaged u-velocity
%           CVavg - Conditionally averaged v-velocity
%           C - Index numbers of set of images used in conditional average

S = size(V);
N = round(S(3)/10); %10-percent of data set
qq = diff(sign(y(1,:)));
yi = (1:length(y(1,:))-1);
zr = yi(logical(qq));
clear qq yi
if SW
    Q = Um(:,zr);  %grab centerline profile
    CI = max(find( Q > max(Q)*0.9));    %find index for end of potential core
    CI = round(CI*0.85);    %reduces core length to ensure we are in core region
    disp(['Correlating in potential core region x/D < ' num2str(x(CI,1))])
    clear Q
    
    [C,I] = max(Um(CI,:));  %Finds peak velocity at end of potential core
    CL = min(find(Um(CI,:)>C/2));
    CL = I - 2*(I-CL);      %Sets top of potential core as center - 2*HWHM
    CR = 2*I -CL;           %Sets bottom of potential core as center + 2*HWHM
else
    CI = length(x(:,1));
    CL = 1;
    CR = S(2);
end
    %crops velocity maps if SW=1 so that correlation only looks at jet potential core
U = U(1:CI,CL:CR,:);
Um = Um(1:CI,CL:CR);
V = V(1:CI,CL:CR,:);
Vm = Vm(1:CI,CL:CR);

q = [];
QV = [];
for n = 1:S(3)
    progress(n,1,S(3),10);
    q = V(:,:,n)-Vm;
    QV(:,n) = q(:);
end
VCM = corrcoef(QV);
clear QV
for n = 1:S(3)
    VCM(n,n) = NaN;
end
done = false;
T = 0.9;
while ~done
    L = sum(VCM > T);
    if max(L) > N
        done = true;
        CorrCoeff = T;
        disp([num2str(max(L)) ' Positive Match Found at Corr. Coeff. = ' num2str(T)])
    elseif T < 0.2
        done = true;
        CorrCoeff = 0;
        disp('No Positive Result')
    else
        T = T-0.01;
    end
end
[M,I] = max(L);
C = unique([I find(VCM(I,:) > T)]);
CVavg = mean(V(:,:,C),3);
CUavg = mean(U(:,:,C),3);

