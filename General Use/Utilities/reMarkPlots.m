function reMarkPlots(h,num)

% reMarkPlots takes a vector of line handles in h
% and reduces the number of plot markers on the lines
% to num. This is useful for closely sampled data.
%
% example:
% t = 0:0.01:pi;
% p = plot(t,sin(t),'-*',t,cos(t),'r-o');
% reMarkPlots(p,10);
% legend('sin(t)','cos(t)')
%
% example 2:
% h = get(gca,'Children');
% reMarkPlots(h,num)
%

for n = 1:length(h)
    if strcmp(get(h(n),'type'),'line')
        axes(get(h(n),'parent'));
        x = get(h(n),'xdata');
        y = get(h(n),'ydata');
        t = 1:length(x);
        s = [0 cumsum(sqrt(diff(x).^2+diff(y).^2))];
        si = (0:num-1)*s(end)/(num-1);
        ti = round(interp1(s,t,si,'linear','extrap'));
        ti = ti+round((rand(1)*2-1)*length(x)/(num-1)/4);
        ti(ti<1) = 1; ti(ti > length(x)) = length(x);
        xi = x(ti);
        yi = y(ti);
        marker = get(h(n),'marker');
        MS = get(h(n),'markersize');
        MEC = get(h(n),'markeredgecolor');
        MFC = get(h(n),'markerfacecolor');
        color = get(h(n),'color');
        style = get(h(n),'linestyle');
        width = get(h(n),'linewidth');
        % make a line with just the markers
        set(line(xi,yi),'marker',marker,'markersize',MS,'markeredgecolor',MEC,...
            'markerfacecolor',MFC,'linestyle','none','linewidth',width,'color',color);
        % make a copy of the old line with no markers
        set(line(x,y),'marker','none','linestyle',style,'linewidth',width,'color',color);
        % set the x- and ydata of the old line to [], this tricks legend to
        % keep on working
        set(h(n),'xdata',xi(1),'ydata',yi(1),'Visible','off');
    end
end