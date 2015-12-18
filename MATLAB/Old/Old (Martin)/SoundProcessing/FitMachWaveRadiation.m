% Assumes the existence of the following variables: 
%   BL_avg - matrix containing baseline narrowband spectra
%   m0 - one dimensional cell array containing matrices of narrowband
%       spectra for each forcing strouhal number. Note that all matrices in 
%       all cases (i.e. m0, m1, etc.) should be the same size. 
%   m1 - same as m0 for azimuthal mode 1
%   m3 - same as m0 for azimuthal mode 3
%   m4 - same as m0 for azimuthal mode \pm4
%   Stdf - a vector containing the forcing strouhal numbers. This of course
%      assumes all modes have the same set of forcing strouhal numbers.


nn = 8; %Channel number
LL = 0.5;   %lower strouhal number bound
UL = 2;     %upper strouhal number bound

%% FIT THE BASELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1;
x = BL_avg(:,2);
y = BL_avg(:,nn);

I = find(x >= LL,1,'first');
Ed = find(x >= UL,1,'first');
if isempty(Ed)
    Ed = length(x);
end

Pval = mmpolyfit(log10(x(I:Ed)),y(I:Ed),1,'Weight',1./x(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
PS(N) = Pval(1);
    %Calculates uncertainty of fit
Serr = sum((Pval(1)*log10(x(I:Ed))+Pval(2)-y(I:Ed)).^2);
merr = length(x(I:Ed));
Derr = merr*sum(log10(x(I:Ed)).^2) - sum(log10(x(I:Ed)))^2;
PE(N) = sqrt(Serr*merr/(merr-2)/Derr);
clear Serr merr Derr Pval x y I Ed
Slope.BL = PS;
SlopeError.BL = PE;
clear PS PE

%% FIT THE FORCED CASES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for N = 1:length(m0)    %N from 1 to number of frequencies
    x = m0{N}(:,2);
    y = m0{N}(:,nn);

    I = find(x >= LL,1,'first');
    Ed = find(x >= UL,1,'first');
    if isempty(Ed)
        Ed = length(x);
    end

    Pval = mmpolyfit(log10(x(I:Ed)),y(I:Ed),1,'Weight',1./x(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
    PS(N) = Pval(1);
        %Calculates uncertainty of fit
    Serr = sum((Pval(1)*log10(x(I:Ed))+Pval(2)-y(I:Ed)).^2);
    merr = length(x(I:Ed));
    Derr = merr*sum(log10(x(I:Ed)).^2) - sum(log10(x(I:Ed)))^2;
    PE(N) = sqrt(Serr*merr/(merr-2)/Derr);
    clear Serr merr Derr Pval x y I Ed
end
Slope.m0 = PS;
SlopeError.m0 = PE;
clear PS PE N

% for N = 1:length(m1)
%     x = m1{N}(:,2);
%     y = m1{N}(:,nn);
% 
%     I = find(x >= LL,1,'first');
%     Ed = find(x >= UL,1,'first');
%     if isempty(Ed)
%         Ed = length(x);
%     end
% 
%     Pval = mmpolyfit(log10(x(I:Ed)),y(I:Ed),1,'Weight',1./x(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
%     PS(N) = Pval(1);
%         %Calculates uncertainty of fit
%     Serr = sum((Pval(1)*log10(x(I:Ed))+Pval(2)-y(I:Ed)).^2);
%     merr = length(x(I:Ed));
%     Derr = merr*sum(log10(x(I:Ed)).^2) - sum(log10(x(I:Ed)))^2;
%     PE(N) = sqrt(Serr*merr/(merr-2)/Derr);
%     clear Serr merr Derr Pval x y I Ed
% end
% Slope.m1 = PS;
% SlopeError.m1 = PE;
% clear PS PE N

for N = 1:length(m3)
    x = m3{N}(:,2);
    y = m3{N}(:,nn);

    I = find(x >= LL,1,'first');
    Ed = find(x >= UL,1,'first');
    if isempty(Ed)
        Ed = length(x);
    end

    Pval = mmpolyfit(log10(x(I:Ed)),y(I:Ed),1,'Weight',1./x(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
    PS(N) = Pval(1);
        %Calculates uncertainty of fit
    Serr = sum((Pval(1)*log10(x(I:Ed))+Pval(2)-y(I:Ed)).^2);
    merr = length(x(I:Ed));
    Derr = merr*sum(log10(x(I:Ed)).^2) - sum(log10(x(I:Ed)))^2;
    PE(N) = sqrt(Serr*merr/(merr-2)/Derr);
    clear Serr merr Derr Pval x y I Ed
end
Slope.m3 = PS;
SlopeError.m3 = PE;
clear PS PE N

for N = 1:length(m4)
    x = m4{N}(:,2);
    y = m4{N}(:,nn);

    I = find(x >= LL,1,'first');
    Ed = find(x >= UL,1,'first');
    if isempty(Ed)
        Ed = length(x);
    end

    Pval = mmpolyfit(log10(x(I:Ed)),y(I:Ed),1,'Weight',1./x(I:Ed)); %Performs a weighted linear fit to account for increased point density in higher frequency
    PS(N) = Pval(1);
        %Calculates uncertainty of fit
    Serr = sum((Pval(1)*log10(x(I:Ed))+Pval(2)-y(I:Ed)).^2);
    merr = length(x(I:Ed));
    Derr = merr*sum(log10(x(I:Ed)).^2) - sum(log10(x(I:Ed)))^2;
    PE(N) = sqrt(Serr*merr/(merr-2)/Derr);
    clear Serr merr Derr Pval x y I Ed
end
Slope.m4 = PS;
SlopeError.m4 = PE;
clear PS PE N

figure
plot(Stdf,Slope.m0,'-sb','LineWidth',1.3)
hold on
% plot(Stdf,Slope.m1,'-or','LineWidth',1.3)
plot(Stdf,Slope.m3,'-d','Color',[102 204 51]/255,'LineWidth',1.3)
plot(Stdf,Slope.m4,'-^k','LineWidth',1.3)
plot(Stdf([1 length(m0)]),[1 1]*Slope.BL,'--k','LineWidth',2.5)
axis([0 3 -35 0])
grid on
xlabel('St_{DF}')
ylabel('(St_{D})^{N}')
title('Semilog Slopes for T_o/T_a = 2.5')
legend('m = 0','m = 3','m = \pm4','BL')
saveas(gcf,'TTR2.5.fig')
saveFigure_v2(gcf,'TTR2.5',600)
% close

save TTR2.5.mat