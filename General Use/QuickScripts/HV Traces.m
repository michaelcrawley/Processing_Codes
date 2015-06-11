n = 7;
TX = '00014';
F = [1 3 5 10 15 20 30];

wave=ReadLeCroyBinaryWaveform(['F:\Research\Jeff Barton Project\20090915-voltage traces\C11' TX '.trc']);
wave.y = nanmoving_average(wave.y,10);
V{n} = [wave.x wave.y];

wave=ReadLeCroyBinaryWaveform(['F:\Research\Jeff Barton Project\20090915-voltage traces\C21' TX '.trc']);
wave.y = nanmoving_average(wave.y,10);
C{n} = [wave.x wave.y];

wave=ReadLeCroyBinaryWaveform(['F:\Research\Jeff Barton Project\20090915-voltage traces\C31' TX '.trc']);
wave.y = nanmoving_average(wave.y,10);
T{n} = [wave.x wave.y];


% subplot(3,1,1); plot(T{n}(:,1)*1e6,T{n}(:,2))
% ylabel('Trigger (V)')
% set(gca,'XTickLabel','');
% axis([-20 10 0 6])
% grid on
% title([num2str(F(n)) ' kHz, 16 \mus, CH 2'])
% 
% subplot(3,1,2); plot(V{n}(:,1)*1e6,V{n}(:,2),'r')
% ylabel('Arc (kV)')
% set(gca,'XTickLabel','');
% axis([-20 10 -5 10])
% grid on
% 
% subplot(3,1,3); plot(C{n}(:,1)*1e6,C{n}(:,2)-0.2369,'k')
% ylabel('Arc (A)')
% axis([-20 10 -1 1])
% grid on
% xlabel('Time (\mus)')
% 
% saveFigure_v2(gcf,[num2str(F(n)) 'kHz'],300)
