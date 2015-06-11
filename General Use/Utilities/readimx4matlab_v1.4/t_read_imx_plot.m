clear

drive ='H';
if isdir([drive ':'])==0
    fprintf('\n====================================================\n');
    fprintf('   Drive [] %s [] is not found\n', upper(drive));
    fprintf('====================================================\n');
end
dname =[drive ':\piv\codes\DaVis_Down\readimx4matlab_v1.4\images\'];
file ='demo2.imx';

addpath =[drive ':\piv\codes\DaVis_Down\readimx4matlab_v1.4\'];

A =readimx([dname file]);
[xx yy V] =showimx(A);

figure(1)
imagesc(V');
axis image
