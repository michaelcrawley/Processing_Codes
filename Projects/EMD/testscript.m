
%Constants
L = 1000;
t = linspace(0,1,L);
dt = t(2)-t(1);
f = 10; %Hz
Kmax = 20;
alpha = 0.25;

%Signals
c = sin(2*pi*f*t); %coherent portion
n = 2*rand(1,L)-1; %noise
p = c+n; %total (noisy) signal

%Initialize
IMF = zeros(Kmax,L);
R = zeros(Kmax,1);

s = p;
for k =1:Kmax
    %Get signal statistics
    avg = wmavg(s,11,'hamming');
    thresh = std(s-avg);
    
    %Find Peaks & Troughs
    [~,P_locs] = findpeaks(s-avg,'minpeakheight',alpha*thresh); %positive peaks
    [~,N_locs] = findpeaks(avg-s,'minpeakheight',alpha*thresh); %negative peaks
    
    %Calculate cubic splines
%     e_max = spline(t(P_locs),[(s(P_locs(1))-s(1))/(t(P_locs(1))-t(1)) s(P_locs) (s(end)-s(P_locs(end)))/(t(end)-t(P_locs(end)))],t);
%     e_min = spline(t(N_locs),[(s(N_locs(1))-s(1))/(t(N_locs(1))-t(1)) s(N_locs) (s(end)-s(N_locs(end)))/(t(end)-t(N_locs(end)))],t);
    e_max = spline([t(1) t(P_locs) t(end)],[0 s(P_locs(1)) s(P_locs) s(P_locs(end)) 0],t);
    e_min = spline([t(1) t(N_locs) t(end)],[0 s(N_locs(1)) s(N_locs) s(N_locs(end)) 0],t);
%     e_max = spline( t(P_locs),[0 s(P_locs) 0],t);
%     e_min = spline(t(N_locs),[0  s(N_locs) 0],t);
    
    %Calculate Trend & IMF
    r = (e_max+e_min)/2; %trend
    IMF(k,:) = s - r;
    s = r;
    
    %Calculate Resemblance Criterion
    if k == 1
        un = IMF(1,:);
        R(k) = 0;
    elseif k > 1 && k < Kmax
        unm = un;
        un = sum(IMF(1:k,:),1);
        R(k) = (un*unm')/(sqrt(un*un')*sqrt(unm*unm'));
        if R(k) < R(k-1)
            IMF(k,:) = IMF(k,:)+r;
            IMF = IMF(1:k,:);
            R = R(1:k);
            break;        
        end
    end
    

end