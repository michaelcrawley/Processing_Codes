clear variables
Sf = 780000;    %Sampling Frequency (Hz)
N = 7800;   %Number of samples per unit

F = 500;   %Forcing Frequency (Hz)
PW = 20e-6; %pulse width (s)
lag = 1e-4; %computation lag (s)

BS = 8192;

t = (0:1/Sf:0.1);  %time axis for one unit
x = 5/2*(square(2*pi*F*t(1:N),F*PW*100)+1);  %pulse train for one unit
NP = t(N)*F;
ph = 2*pi*(NP-floor(NP));

% plot(t(1:N),x)

done = false;
q = min(find(t>t(length(x))+1e-5+rand*lag));
while ~done
    x = [x ones(1,q-length(x)-1)*x(end) 5/2*(square(2*pi*F*t(1:N)+ph,F*PW*100)+1)];
    if length(x) > length(t)
        done = true;
        x = x(1:length(t));
    else
        NP = t(length(x))*F;
        ph = 2*pi*(NP-floor(NP));
        q = min(find(t>t(length(x))+lag+rand*1e-4));
    
        if isempty(q)
            done = true;
            x = [x ones(1,length(t)-length(x))*x(end)];
        end
    end
        
%     plot(t(1:length(x)),x)
end
% plot(t,x)
    
daq = interp1(t,x,(0:1/200000:0.1));

figure
for n = 1:100
    S = round(rand(1)*(length(daq)-BS-1))+1;
    [f,mx(:,n),px]=myfft(daq(S:S+BS-1),200000);
%     m = sqrt(mean(mx.^2,2));
%     loglog(f,m)
%     axis([100 100000 1e-4 1])
%     grid on
%     pause(0.1)
end
m = sqrt(mean(mx.^2,2));
loglog(f,m)
axis([100 100000 1e-4 1])
grid on    


