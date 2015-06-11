files = who;

az = '90';

%Calc ICAO values
Tref = 25+273;
Te =  0.68965517*Tref;
Uref = 1.5*sqrt(1.4*287*Te);

%Pull ambient conditions
ChPol =eval([files{1},'.processing_params.ChPol']);
f = eval([files{1},'.processing_params.fx']);
df = mean(diff(f));
D = eval([files{1},'.processing_params.D']);
Ta =  eval([files{1},'.processing_params.AT'])+273;
a = sqrt(1.4*287*Ta);

for n = 1:length(files),
    U(n) = eval([files{n},'.data.Ue']);
    Mj(n) = U(n)/a; Mc(n) = 0.65*Mj(n);
    SPL(:,:,n) = eval([files{n},'.data.dBCorrected']);
    for nn = 1:length(ChPol)
       SPL_s(:,nn,n) = SPL(:,nn,n)-10*log10(df)-10*log10(D./U(n))-80*log10(U(n)/Uref)+50*log10(1-Mc(n)*cosd(ChPol(nn)));
       st(:,nn,n) = (f*D*(1-Mc(n)*cosd(ChPol(nn)))/U(n));
    end    
end

%calc averages
avgSPL = mean(SPL,3);
avgSPL_s = mean(SPL_s,3);

save(['PSD data az' az],'Tref','Uref','U','SPL','files','ChPol','D','SPL_s','avgSPL','avgSPL_s','st','f');