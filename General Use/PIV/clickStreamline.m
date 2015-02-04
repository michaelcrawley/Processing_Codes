function clickStreamline(x,y,U,V,xEsc,color)
%clickStreamline(x,y,U,V,xEsc,color) 

hold on
[qx,qy] = ginput(1);
while qx < xEsc
    h = streamline(x',y',U',V',qx,qy);
    set(h,'Color',color)
    [qx,qy] = ginput(1);
end
hold off
