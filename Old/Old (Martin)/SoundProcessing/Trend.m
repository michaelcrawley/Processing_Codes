function [slope,slopeError,T] = Trend(N,ChN);

% N = 21;  %Number of Spectra
% ChN = 6;  %Data Column from which to generate slopes

colors = 'krbgmc';
style{1} = '-';
style{2} = '--';
style{3} = '-.';
style{4} = ':';
done2 = false;
toto = 0;
n = 0;
while ~done2
    [SPfile,SPpath] = uigetfile('*.fftNOS',['Select Spectrum: ' num2str(n+1)],'MultiSelect','on');
    if ~iscell(SPfile)
        stemp = SPfile;
        clear SPfile
        SPfile{1} = stemp;
        clear stemp
    end
    totn = toto+length(SPfile);
    for n = toto+1:totn
        Dat{n}=dlmread([SPpath SPfile{n-toto}],'\t',1,0);
        if exist([SPpath SPfile{n-toto}(1:end-6) 'SET'],'file')
            done = false;
            fid = fopen([SPpath SPfile{n-toto}(1:end-6) 'SET'],'r');
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Jet Temperature Ratio:',L));
                if mtch
                    T{n} = L(24:end);
                    done = true;
                end
            end
            fclose(fid);
        elseif exist([SPpath SPfile{n-toto}(1:end-9) 'SET'],'file')
            done = false;
            fid = fopen([SPpath SPfile{n-toto}(1:end-9) 'SET'],'r');
            while ~done
                L = fgetl(fid);
                mtch = ~isempty(strmatch('Jet Temperature Ratio:',L));
                if mtch
                    T{n} = L(24:end);
                    done = true;
                end
            end
            fclose(fid);
        else
            T{n} = input(['  File ' num2str(n) ': Enter Temp Ratio: '],'s');
        end
        Dat{n} = Dat{n}(10:end-5,:);
        semilogx(Dat{n}(:,2),Dat{n}(:,ChN),[style{mod(ceil(n/6),4)+1} colors(mod(n,6)+1)])
        hold on
    end
    if totn >= N
        done2 = true;
        N = totn;
    else
        toto = totn;
    end
end
grid on
legend(T)

slope = [];
for n = 1:N
    T{n} = str2num(T{n});
    data = Dat{n};
    [c,I] = max(data(:,ChN));
    P = mmpolyfit(log(data(I*2:round(end/4),2)),data(I*2:round(end/4),ChN),1,'Weight',1./data(I*2:round(end/4),2));
    S = sum((P(1)*log(data(I*2:round(end/4),2))+P(2)-data(I*2:round(end/4),ChN)).^2);
    m = length(data(I*2:round(end/4),2));
    D = m*sum(log(data(I*2:round(end/4),2)).^2) - sum(log(data(I*2:round(end/4),2)))^2;
    PE = sqrt(S*m/(m-2)/D);
    
    semilogx(data(I:end,2),P(1)*log(data(I:end,2))+P(2),'--k')
    text(data(I+10,2),P(1)*log(data(I+10,2))+P(2),['St_D^{' num2str(P(1),'%10.3f') ' \pm ' num2str(PE,'%10.3f') '}'])
    slope(n) = P(1);
    slopeError(n) = PE;
end
hold off
xlabel('St_D')
ylabel('SPL')
T = cell2mat(T);