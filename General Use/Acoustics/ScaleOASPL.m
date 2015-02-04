files = who;

az = 'M130';

%Calc ICAO values
TICAO = 25+273;
PICAO = 101325;
rhoICAO = PICAO/287.058/TICAO;
Te =  0.68965517*(25+273);
UICAO = 1.5*sqrt(1.4*287*Te);

%Pull .mat file data
for n = 1:length(files),
    rho(n) = eval([files{n},'.data.rhoe']);
    U(n) = eval([files{n},'.data.Ue']); 
    OASPL(n,:) = eval([files{n},'.results.OASPL']);
end

%Pull ambient conditions
ChPol =eval([files{1},'.processing_params.ChPol']);
Ta = eval([files{1},'.processing_params.AT'])+273;
Pa = eval([files{1},'.processing_params.AP']);
rhoa = Pa/287.058/Ta*1000;

%calc scaled OASPL
for n = 1:length(files)
    OASPL_s(n,:) = OASPL(n,:)-10*log10(((U(n)/UICAO)^8)*(rho(n)*rhoa/rhoICAO/rhoICAO)*(Ta/TICAO)^2);
end

%calc averages
avgOASPL = mean(OASPL,1);
avgOASPL_s = mean(OASPL_s,1);

h = figure;
plot(ChPol,OASPL); 
xlabel('Polar Angle');xlim([25 130]);
ylabel('OASPL (dB)');ylim([108 122]);
set(gca,'XDir','reverse');
legend(files,'Interpreter','none','Location','Northwest');
title(['Unscaled OASPL for Az: ' az '^o']);
grid on;
saveas(h,['OASPL az' az],'fig');
saveas(h,['OASPL az' az],'png');
close(h);

h = figure;
plot(ChPol,OASPL_s); 
xlabel('Polar Angle');xlim([25 130]);
ylabel('Scaled OASPL (dB)');ylim([108 122]);
set(gca,'XDir','reverse');
legend(files,'Interpreter','none','Location','Northwest');
title(['Scaled OASPL for Az: ' az '^o']);
grid on;
saveas(h,['scaled OASPL az' az],'fig');
saveas(h,['scaled OASPL az' az],'png');
close(h);

save(['OASPL data az' az],'TICAO','PICAO','rhoICAO','UICAO','rho','U','OASPL','files','ChPol','Ta','Pa','OASPL_s','avgOASPL','avgOASPL_s');