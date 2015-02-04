function [Matches,CorrCoeff,SS] = PIV_CorrMap_v3(x,y,Um,V,Vm,SW)
% % % function [CUavg,CVavg,Matches,CorrCoeff,SS] = PIV_CorrMap_v3(x,y,U,Um,V,Vm,SW)
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
%               - else == Correlate over this distance
%
%Outputs:   Matches - Index numbers of set of images used in conditional average
%           CorrCoeff - Correlation coefficent of matches - zero if no match
%           SS - Structure spacing - zero if no match

LL = 0.2;   %Acceptable limit on correlation coefficient


S = size(V);
N = round(S(3)/10); %10-percent of data set
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
        CO = min(SW); CI = max(SW);
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
    %crops velocity maps if SW=1 so that correlation only looks at jet potential core
% U = U(CO:CI,CL:CR,:);
% Um = Um(CO:CI,CL:CR);
V = V(CO:CI,CL:CR,:);
Vm = Vm(CO:CI,CL:CR);

q = []; QV = [];
for n = 1:S(3)  %Processes velocity fields into form for "corrcoef" function
    progress(n,1,S(3),10);
    q = V(:,:,n)-Vm;    %Removes mean velocity
    QV(:,n) = q(:);
end
VCM = corrcoef(QV);
% clear QV
for n = 1:S(3)  %Removes primary diagonal elements
    VCM(n,n) = NaN;
end

done = false; threshold = 0.9;
while ~done
    L = sum(VCM > threshold);   %Finds the number of images with correlation above "threshold"
    if max(L) > N   %If the necessary number of images is found
        done = true;
        CorrCoeff = threshold;
        disp([num2str(max(L)) ' Positive Match Found at Corr. Coeff. = ' num2str(threshold)])
    elseif threshold < LL   %If threshold has fallen below acceptable lower limit
        done = true;
        CorrCoeff = 0;
        disp('No Positive Result')
    else
        threshold = threshold-0.01;
    end
end
[M,I] = max(L); %"I" is index of image used as reference for first phase
Matches = unique([I find(VCM(I,:) > threshold)]); %Store indices of images which meet condition
CVavg = mean(V(:,:,Matches),3);
% CUavg = mean(U(:,:,Matches),3);


done = false; threshold = -0.9;
while ~done
    L = sum(VCM(I,:) < threshold);  %Looks for negative correlation match to reference image "I"
    if L > N
        done = true;
        CorrCoeff(2) = threshold;
        disp([num2str(L) ' Negative Match Found at Corr. Coeff. = ' num2str(threshold)])
    elseif threshold > -LL
        done = true;
        CorrCoeff(2) = 0;
        disp('No Negative Result')
    else
        threshold = threshold+0.01;
    end
end
M = find(VCM(I,:) < threshold);
Matches(2,1:length(M)) = M;   %Store indices of images with negative correlation
CVavg(:,:,2) = mean(V(:,:,M),3);
% CUavg(:,:,2) = mean(U(:,:,M),3);


XC = xcorr2(CVavg(:,:,1),CVavg(:,:,2)); %Compute correlation of two phases - since two phases are pi separated, max negative correlation should occur at zero.
I = ceil(size(XC,2)/2); XC = XC(:,I);   %Isolate the y' = 0 correlation
I = ceil(size(XC,1)/2);
% xC = ((1:length(XC))'-I); %Create grid coordinates for correlation
dXC = [0; diff(XC)];
NP = find(dXC(I+3:end)<0,1,'first')+I+1 -find(dXC(1:I-3)>0,1,'last');
%     I = ezfit(xC,XC/max(XC),'-cos(a*x).*exp(-x.^2/b)'); %Fit the correlation function with a Gaussian windowed cosine function
%     NP = abs(2*pi/I.m(1));   %Find number of points in structure period.

    %Compute 3rd and 4th phases if first two phases produced results
if sum(abs(CorrCoeff))>0
%     CUavg(:,:,3) = CUavg(:,:,2); CUavg(:,:,2) = 0;
    CVavg(:,:,3) = CVavg(:,:,2); CVavg(:,:,2) = 0;
    Matches(3,:) = Matches(2,:); Matches(2,:) = 0;
    CorrCoeff(3) = CorrCoeff(2); CorrCoeff(2) = 0;

    L2 = round(NP/4)+1;
    CV2 = zeros(size(CVavg(:,:,1)));
    CV2(L2:end,:) = CVavg(1:end-L2+1,:,1)-Vm(1:end-L2+1,:); %Create flow field with structures shifted by L/4

    CORR = zeros(S(3),1); MX = sum(CV2(:).^2);
    for n = 1:length(CORR)  %Correlate images with shifted structure image
        CORR(n) = sum(CV2(:).*QV(:,n))/MX;
    end

    done = false;
    threshold = 1;
    while ~done
        L = sum(CORR > threshold);
        if L > N
            done = true;
            CorrCoeff(2) = threshold;
            disp([num2str(max(L)) ' 3rd Phase Match Found at Corr. Coeff. = ' num2str(threshold)])
        elseif threshold < LL
            done = true;
            CorrCoeff(2) = 0;
            disp('No 3rd Phase Result')
        else
            threshold = threshold-0.01;
        end
    end
    M = find(CORR > threshold);
    Matches(2,1:length(M)) = M;
%     CVavg(:,:,2) = mean(V(:,:,M),3);
%     CUavg(:,:,2) = mean(U(:,:,M),3);


    done = false;
    threshold = -1;
    while ~done
        L = sum(CORR < threshold);
        if L > N
            done = true;
            CorrCoeff(4) = threshold;
            disp([num2str(L) ' 4th Phase Match Found at Corr. Coeff. = ' num2str(threshold)])
        elseif threshold > -0.2
            done = true;
            CorrCoeff(4) = 0;
            disp('No Result')
        else
            threshold = threshold+0.01;
        end
    end
    M = find(CORR < threshold);
    Matches(4,1:length(M)) = M;
%     CVavg(:,:,4) = mean(V(:,:,M),3);
%     CUavg(:,:,4) = mean(U(:,:,M),3);
end

SS = NP*mean(diff(x(:,1)));  %Structure spacing - units of x