function [LS,fL,A] = LSS(f,SPL)

% [FileName,PathName] = uigetfile('*.fftNOS','Select Spectrum File','MultiSelect','off');
% data = dlmread([PathName FileName],'\t', 1, 0);

% fC = input('  Enter Frequency Column Number: ');
% dC = input('  Enter Data Column Number: ');
% f = data(:,fC);
% SPL = data(:,dC);

[A,I] = max(SPL);   %location of spectrum peak
fL = f(I);
done = false;

figure;

while ~done
    fofL = f/fL;

    L = zeros(size(SPL));

    L(fofL <= 0.5) = 2.53895+18.4*log10(fofL(fofL <= 0.5));

    L((fofL <= 1)&(fofL > 0.5)) = -38.19338*log10(fofL((fofL <= 1)&(fofL > 0.5))).^2-...
        16.91175*log10(fofL((fofL <= 1)&(fofL > 0.5))).^3;

    L((fofL <= 2.5)&(fofL > 1)) = (1.06617 - 45.2994*log10(fofL((fofL <= 2.5)&(fofL > 1)))+...
        21.40972*log10(fofL((fofL <= 2.5)&(fofL > 1))).^2).*log10(fofL((fofL <= 2.5)&(fofL > 1)));

    L(fofL > 2.5) = 5.64174-27.7472*log10(fofL(fofL > 2.5));

    LS = A+L;

    semilogx(f,SPL)
    hold on
    semilogx(f,LS,'r')
    hold off
    grid on
    
    fprintf(['  A: ' num2str(A) ' \n'])
    fprintf(['  F: ' num2str(fL) '\n'])
    ok = input('  Is this Good (y/n): ','s');
    if strcmp(ok,'y')
        done = true;
    else
        A = input('  Enter Amplitude: ');
        fL = input('  Enter Peak Frequency: ');
    end
end
title({'Large Scale Similiarity Spectrum', ['A = ' num2str(A) '  fL = ' num2str(fL)]})

% [pathstr,name,ext] = fileparts(FileName);
% saveas(gcf,[PathName name '_LSS'],'fig');
% saveas(gcf,[PathName name '_LSS'],'png');
% 
% fprintf('Done\n');
