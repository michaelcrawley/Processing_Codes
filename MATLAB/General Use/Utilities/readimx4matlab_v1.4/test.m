clear



drive ='H';
if isdir([drive ':'])==0
    fprintf('\n====================================================\n');
    fprintf('   Drive [] %s [] is not found\n', upper(drive));
    fprintf('====================================================\n');
end
dname =[drive ':\jhkim\piv\20060606\T09_0030_74_05_Vm11\Vel_32_50\'];
file ='B00001.vc7';



addpath =[drive '\codes\DaVis_Down\readimx4matlab_v1.4\'];


for kf=1:10000
    file =sprintf('B%05d.vc7',kf);
    tmp =fopen([dname file],'r');
    if tmp ==-1
        break;
    else
        fclose(tmp);
    end
end
num_files =kf-1;

for kf=1:num_files
    if mod(kf-1,10)==0
        skip =num2str(kf);
        fprintf('%s',skip(end));
    end
    file =sprintf('B%05d.vc7',kf);
    tmp =fopen([dname file],'r');
    if tmp ==-1
        break;
    else
        fclose(tmp);
    end
    A =readimx([dname file]);
    [xx yy V_x V_y] =showimx(A);
    if kf==1
        sum_x =V_x;
        sum_y =V_y;
    else
        sum_x =sum_x+V_x;
        sum_y =sum_y+V_y;
    end
end
fprintf('\n');
num_files =kf-1;
avg_V_x =sum_x/num_files;
avg_V_y =sum_y/num_files;

[num_xx num_yy]=size(avg_V_x);

figure(1)
imagesc(avg_V_x');

figure(2)
contourf(xx,yy,avg_V_x,10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Save in DAT format
%%
%%%%%%%%%%%%%%%
file_dat =[dname 'VEL_AVG.DAT'];
dat =fopen(file_dat,'w');

if dname(end)~='\'
    dname =[dname '\'];
end
id_slash =strfind(dname,'\');
title_dat =dname(id_slash(end-2)+1:id_slash(end-1)-1);
    
fprintf(dat,'TITLE = "%s"\n',title_dat);
fprintf(dat,'VARIABLES = "x", "y", "Vx", "Vy"\n');
fprintf(dat,'ZONE I=%g, J=%g, F=POINT\n',num_xx,num_yy);
for kp =1:num_yy
    for jp =1:num_xx
        fprintf(dat,'%g %g %g %g\n',xx(jp,kp),yy(jp,kp),avg_V_x(jp,kp),avg_V_y(jp,kp));
    end
end
fclose(dat);
