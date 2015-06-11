function [IMF R] = EMD(p)


    %Constants
    Kmax = 20;
    alpha = 0.125;
    L = length(p);

    %Initialize
    p = p(:);
    IMF = zeros(L,Kmax);
    R = zeros(Kmax,1);

    s = p;
    for k =1:Kmax
        %Get signal statistics
        avg = wmavg(s,5,'rectwin',1);
        thresh = std(s-avg);

        %Find Peaks & Troughs
        [~,P_locs] = findpeaks(s-avg,'minpeakheight',alpha*thresh); %positive peaks
        [~,N_locs] = findpeaks(avg-s,'minpeakheight',alpha*thresh); %negative peaks

        %Calculate cubic splines
        e_max = spline([1; P_locs; L],[0; s(P_locs(1)); s(P_locs); s(P_locs(end)); 0],1:L);
        e_min = spline([1; N_locs; L],[0; s(N_locs(1)); s(N_locs); s(N_locs(end)); 0],1:L);

        %Calculate Trend & IMF
        r = (e_max(:)+e_min(:))/2; %trend
        IMF(:,k) = s - r;
        s = r;

        %Calculate Resemblance Criterion
        if k == 1
            un = IMF(:,1);
            R(k) = 0;
        elseif k > 1 && k < Kmax
            unm = un;
            un = sum(IMF(:,1:k),2);
            R(k) = (un'*unm)/(sqrt(un'*un)*sqrt(unm'*unm));
            if R(k) < R(k-1)
                IMF(:,k) = IMF(:,k)+r;
                IMF = IMF(:,1:k);
                R = R(1:k);
                break;        
            end
        end
    end

end