function yTickSeparationMarker(dx,dy,zx,zy)
%This function prints a yTick spacing marker for plots in which it is
%desired to have no absolute scale. The output includes a graphical object
%similar to "]" and the height of said bar (e.g. "]5")
%
%Function Call: yTickSeparationMarker(dx,dy,zx,zy)
%
%INPUTS:
%  dx - the abscissa extent of the marker
%  dy - the ordinate extent of the marker
%  zx - the abscissa coordinate of the lower left corner of the marker
%  zy - the ordinate coordinate of the lower left corner of the marker


C = [0 0 0]/255;

hold on
plot([zx zx+dx],[zy zy],'LineWidth',2,'Color',C);
plot([zx zx+dx],[zy zy]+dy,'LineWidth',2,'Color',C);
plot([zx zx]+dx,[zy zy+dy],'LineWidth',2,'Color',C);
text(zx+dx,zy+dy/2,[' ' num2str(dy)],'Color',C);

hold off
