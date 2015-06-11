function minorGridlineFix(h)
% function minorGridlineFix(h,XP,YP,LS)

AP = get(h);
h2 = axes;

% % AP2 = AP;
% % AP2 = rmfield(AP2,{'BeingDeleted' ;'Children' ; });



set(h2,'Color','none')
set(h2,'XTickLabel','')
set(h2,'YTickLabel','')
set(h2,'Position',AP.Position)
set(h2,'XDir',AP.XDir)
set(h2,'YScale',AP.YScale)
set(h2,'XLim',AP.XLim)
set(h2,'YLim',AP.YLim)

set(h2,'GridLineStyle','--');
set(h2,'MinorGridLineStyle','--');
axes(h2); grid on


% % % %
% % % %
% % % %FUNCTION CALLS
% % % % minorGridlineFix(gca,XP,[],LS)
% % % % minorGridlineFix(gca,[],YP,LS)
% % % % minorGridlineFix(gca,XP,YP,LS)
% % % 
% % % 
% % % AP = get(h);
% % % 
% % % if ~isempty(XP)
% % %     set(h,'XMinorGrid','off')
% % %     set(h,'XMinorTick','off')
% % %     
% % %     T = AP.XTick; 
% % %     L = AP.XLim;
% % %     if strcmp(AP.XScale,'linear')
% % %         dT = T(2)-T(1); LR = [floor(L(1)/dT) ceil(L(2)/dT)]*dT;
% % %         MG = (LR(1):dT/XP:LR(2));
% % %     else
% % %         OT = log10(T); OL = [floor(log10(L(1))) ceil(log10(L(2)))];
% % %         OM = 10.^(union(OT,OL));
% % %         for n = 1:length(OM)-1
% % %             MG = [MG linspace(OM(n),OM(n+1),XP)];
% % %         end
% % %         MG = MG(MG > L(1)); MG = MG(MG < L(2));
% % %     end
% % %     MG = unique(MG);
% % %     MG = setdiff(MG,T);
% % %     
% % %     hold on
% % %     for n = 1:length(MG)
% % %         plot(h,[MG(n) MG(n)],AP.YLim,LS)
% % %     end
% % %     hold off
% % % end
% % % 
% % % if ~isempty(YP)
% % %     set(h,'YMinorGrid','off')
% % %     set(h,'YMinorTick','off')
% % %     
% % %     T = AP.YTick;
% % %     L = AP.YLim;
% % %     if strcmp(AP.YScale,'linear')
% % %         dT = T(2)-T(1); LR = [floor(L(1)/dT) ceil(L(2)/dT)]*dT;
% % %         MG = (LR(1):dT/YP:LR(2));
% % %     else
% % %         OT = log10(T); OL = [floor(log10(L(1))) ceil(log10(L(2)))];
% % %         OM = 10.^(union(OT,OL));
% % %         MG = [];
% % %         for n = 1:length(OM)-1
% % %             MG = [MG linspace(OM(n),OM(n+1),YP+2)]; (OM(n):(OM(n+1)-OM(n))/YP:OM(n+1))
% % %         end
% % %         MG = MG(MG > L(1)); MG = MG(MG < L(2));
% % %     end
% % %     MG = unique(MG);
% % %     MG = setdiff(MG,T);
% % %     
% % %     hold on
% % %     for n = 1:length(MG)
% % %         plot(h,AP.XLim,[MG(n) MG(n)],LS)
% % %     end
% % %     hold off
% % % end
