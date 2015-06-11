function [FS,fF,B] = FSS(f,SPL)

% [FileName,PathName] = uigetfile('*.fftNOS','Select Spectrum File','MultiSelect','off');
% data = dlmread([PathName FileName],'\t', 1, 0);
% 
% fC = input('  Enter Frequency Column Number: ');
% dC = input('  Enter Data Column Number: ');
% f = data(11:end,fC);
% SPL = data(11:end,dC);

[B,I] = max(SPL);   %location of spectrum peak
fF = f(I);
done = false;

figure;

while ~done
    fofF = f/fF;

    F = zeros(size(SPL));

    F(fofF <= 0.05) = 9.9+14.91126*log10(fofF(fofF <= 0.05));

    F((fofF <= 0.15)&(fofF > 0.05)) = -3.5 + (11.874876 + 2.1202444*log10(20/3*fofF((fofF <= 0.15)&(fofF > 0.05)))+...
        7.5211814*log10(20/3*fofF((fofF <= 0.15)&(fofF > 0.05))).^2).*log10(20/3*fofF((fofF <= 0.15)&(fofF > 0.05)));

    F((fofF <= 1)&(fofF > 0.15)) = (-1.0550362 + 4.9774046*log10(fofF((fofF <= 1)&(fofF > 0.15)))).*log10(fofF((fofF <= 1)&(fofF > 0.15))).^2;

    F((fofF <= 10)&(fofF > 1)) = -(8.1476823 + 3.6523177*log10(fofF((fofF <= 10)&(fofF > 1)))).*log10(fofF((fofF <= 10)&(fofF > 1))).^2;

    F((fofF <= 30)&(fofF > 10)) = -11.8 - (27.2523 + 0.8091863*log10(fofF((fofF <= 30)&(fofF > 10))/10)+...
        14.851964*log10(fofF((fofF <= 30)&(fofF > 10))/10).^2).*log10(fofF((fofF <= 30)&(fofF > 10))/10);

    F(fofF > 30) = 29.77786-38.16739*log10(fofF(fofF > 30));

    FS = B+F;

    semilogx(f,SPL)
    hold on
    semilogx(f,FS,'r')
    hold off
    grid on
    
    fprintf(['  B: ' num2str(B) ' \n'])
    fprintf(['  F: ' num2str(fF) '\n'])
    ok = input('  Is this Good (y/n): ','s');
    if strcmp(ok,'y')
        done = true;
    else
        B = input('  Enter Amplitude: ');
        fF = input('  Enter Peak Frequency: ');
    end
end
title({'Fine Scale Similiarity Spectrum', ['B = ' num2str(B) '  fF = ' num2str(fF)]})

% [pathstr,name,ext] = fileparts(FileName);
% saveas(gcf,[PathName name '_FSS'],'fig');
% saveas(gcf,[PathName name '_FSS'],'png');
% 
% fprintf('Done\n');
