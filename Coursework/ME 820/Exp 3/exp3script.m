clear;clc;

fname = 'experiment3.xls';
sheetname = {'Straight Pipe (1000 RPM, WOT)', 'Straight Pipe (2000 RPM, WOT)', 'Muffler (1000 RPM, WOT)','Muffler (2000 RPM, WOT)'};
tmp = zeros(7200,3,4);
RPM = [1000 2000 1000 2000];
for j = 1:length(sheetname)
   tmp(:,:,j) = xlsread(fname,sheetname{j});
   time{j} = (tmp(:,1,j)/360/RPM(j))*60;%check this
   FS = mean(diff(time{j}))^-1;
   N = length(time{j});
   [f{j}, SPL{j}] = calcSPL_v2(tmp(:,3,j),FS,N);
   SPL{j} = 10*log10(SPL{j}/4E-20);
   O{j} = f{j}*60/RPM(j);
end

%Part 1
h(1) = figure; j = 3;
plot(time{j},tmp(:,2,j),'k',time{j},tmp(:,3,j),'k--');title(sheetname{j});xlabel('Time (s)');ylabel('Pressure (bar)');legend('Upstream','Downstream');
filename = sheetname{j};
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

h(2) = figure; j = 4;
plot(time{j},tmp(:,2,j),'k',time{j},tmp(:,3,j),'k--');title(sheetname{j});xlabel('Time (s)');ylabel('Pressure (bar)');legend('Upstream','Downstream');
filename = sheetname{j};
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%Part 2
h(3) = figure; 
plot(time{1},tmp(:,3,1),'k',time{3},tmp(:,3,3),'k--');title('1000 RPM Pressure');xlabel('Time (s)');ylabel('Pressure (bar)');legend('Straight Pipe','Muffler');
filename = '1000 RPM Pressure';
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

h(4) = figure; 
plot(time{2},tmp(:,3,2),'k',time{4},tmp(:,3,4),'k--');title('2000 RPM Pressure');xlabel('Time (s)');ylabel('Pressure (bar)');legend('Straight Pipe','Muffler');
filename = '2000 RPM Pressure';
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%Part 3
h(5) = figure;
plot(O{1},SPL{1},'k',O{3},SPL{3},'k--');title('1000 RPM Pressure Spectra');legend('Straight Pipe','Muffler');xlim([0 20]);xlabel('Order');ylabel('SPL (dB)');
filename = '1000 Spectra';
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');
fid = fopen([filename,'.txt'],'w');
fprintf(fid,'%6.1f \t %6.1f \t %6.1f \n',[O{1} SPL{1} SPL{3}]');
fclose(fid);

h(6) = figure;
plot(O{2},SPL{2},'k',O{4},SPL{4},'k--');title('2000 RPM Pressure Spectra');legend('Straight Pipe','Muffler');xlim([0 20]);xlabel('Order');ylabel('SPL (dB)');
filename = '2000 Spectra';
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');
fid = fopen([filename,'.txt'],'w');
fprintf(fid,'%6.1f \t %6.1f \t %6.1f \n',[O{2} SPL{2} SPL{4}]');
fclose(fid);

%Part 4
fid = fopen('Insertion loss.txt','w');
fprintf(fid,'%6.1f \t %6.3f \n',[1000 SPL{1}(5)-SPL{3}(5); 2000 SPL{2}(5)-SPL{4}(5)]');
fclose(fid);

close(h);





