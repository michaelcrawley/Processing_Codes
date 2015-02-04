clear all;
clc;

filename = 'Experiment1.xls';
sheetname = {'0.1'; '0.15'; '0.2' ;'0.25' ;'0.3'};
range = 'A1:D752';
data = cell(1,length(sheetname)+1);
cutoff = 7;
Temp =  25+273; %in K
c = sqrt(287*1.4*Temp); % in m/s
D = 0.6*0.0254; %in m
lc = 0.0254*[1 8]+0.85*D/2;

marker = [ 'x'; '+'; '*'; 'p'; '^'];
fid1 = fopen('l1amplitudes.txt','w');
fid2 = fopen('l8amplitudes.txt','w');

h(1) = figure;title('Frequency of Whistle Tones vs Mach Number for L = 1 in');xlabel('Mach Number');ylabel('Frequency (Hz)');hold on;
h(2) = figure;title('Frequency of Whistle Tones vs Mach Number for L = 8 in');xlabel('Mach Number');ylabel('Frequency (Hz)');hold on;
h(3) = figure;title('Amplitude of Whistle Tones vs Strouhal Number for L = 1 in');xlabel('Strouhal Number');ylabel('Amplitude (dB)');hold on;
h(4) = figure;title('Amplitude of Whistle Tones vs Strouhal Number for L = 8 in');xlabel('Strouhal Number');ylabel('Amplitude (dB)');hold on;

for i = 1:length(sheetname)
    tmp = xlsread(filename,sheetname{i},range);
    noise.L0 = mean(tmp(:,2));
    noise.L1 = mean(tmp(:,3));
    noise.L8 = mean(tmp(:,4));
    
    [~, loc.L0] = findpeaks(tmp(:,2),'MINPEAKHEIGHT',noise.L0+cutoff);
    [~, L1] = findpeaks(tmp(:,3),'MINPEAKHEIGHT',noise.L1+cutoff);
    [~, L8] = findpeaks(tmp(:,4),'MINPEAKHEIGHT',noise.L8+cutoff);
    loc.L1 = [];
    loc.L8 = [];
    
    for n = 1:length(L1)
        test1 = any(abs(L1(n)-loc.L0) < 50); %baseline peak test
        test2 = abs(L1(n)-L1)<50; %neighbor peak test
        if ~test1 && tmp(L1(n),3) == max(tmp(L1(test2),3))
            loc.L1 = [loc.L1 L1(n)];
            fprintf(fid1,'%6.2f \t %6.0f \t %6.2f \t %6.2f \n',[str2double(sheetname{i}) tmp(L1(n),1) tmp(L1(n),1)*D/c/str2double(sheetname{i}) tmp(L1(n),3)]);
        end
    end
    
    for n = 1:length(L8)
        test1 = any(abs(L8(n)-loc.L0) < 50); %baseline peak test
        test2 = abs(L8(n)-L8)<50; %neighbor peak test
        if ~test1 && tmp(L8(n),4) == max(tmp(L8(test2),4))
            loc.L8 = [loc.L8 L8(n)];
            fprintf(fid2,'%6.2f \t %6.0f \t %6.2f \t %6.2f \n',[str2double(sheetname{i}) tmp(L8(n),1) tmp(L8(n),1)*D/c/str2double(sheetname{i}) tmp(L8(n),4)]);
        end
    end
        
    figure(h(1));
    plot(str2double(sheetname{i})*ones(1,length(loc.L1)),tmp(loc.L1,1),'*');
    figure(h(2));
    plot(str2double(sheetname{i})*ones(1,length(loc.L8)),tmp(loc.L8,1),'*');
    figure(h(3));
    plot(tmp(loc.L1,1)*D/c/str2double(sheetname{i}), tmp(loc.L1,3),marker(i));
    figure(h(4));
    plot(tmp(loc.L8,1)*D/c/str2double(sheetname{i}), tmp(loc.L8,4),marker(i));
    
    data{i} =  struct('Mach_number',sheetname{i},'Frequency',tmp(:,1),'L0in',tmp(:,2),'L1in',tmp(:,3),'L8in',tmp(:,4),'Peaks',loc); 
    sheetname{i} = ['M = ' sheetname{i}];
end

figure(h(3)); legend(sheetname');
figure(h(4)); legend(sheetname');

saveas(h(1),'FvM1','fig');saveas(h(1),'FvM1','png');
saveas(h(2),'FvM8','fig');saveas(h(2),'FvM8','png');
saveas(h(3),'AvS1','fig');saveas(h(3),'AvS1','png');
saveas(h(4),'AvS8','fig');saveas(h(4),'AvS8','png');

fclose(fid1);fclose(fid2);
close(h);

ft1 = (2*(0:6)+1)*c/4/lc(1);
ft2 = (2*(0:6)+1)*c/4/lc(2);
fid = fopen('theoretical frequencies.txt','w');
fprintf(fid,'%6.0f \t %6.0f \n',[ft1; ft2]);
fclose(fid);

tmp = xlsread(filename,'Broadband',range);
data{end} = struct('Mach_number','Broadband','Frequency',tmp(:,1),'Straight',tmp(:,2),'Sharp_Turn',tmp(:,3)); 

h = figure;
plot(tmp(:,1),tmp(:,2),tmp(:,1),tmp(:,3));title('Broadband noise generation'); xlabel('Frequency (Hz)'); ylabel('SPL (dB)');legend('Straight duct','Bent duct');
saveas(gcf,'broadband','fig');saveas(gcf,'broadband','png');
close(h);