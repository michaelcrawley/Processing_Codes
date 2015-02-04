function axisEqual(AX)
%This function equalizes the axes of a graph, similar to the command 
%"axis equal", but does so without the strange sizing artifacts which screw 
%up image extent when printed to file. NOTE: CURRENT VERSION IS FOR 
%TWO-DIMENSIONAL AXES ONLY AND DOES NOT WORK FOR SUBPLOTS. Also, if you 
%intend to use a colorbar or axis labels which take up a lot of space, this 
%command should be called  after the creation of these objects because 
%invoking them will distort the axes.
%
%INPUT
% AX - the handle for the axis to be manipulated

UN = get(AX,'Units');   %Gets the axis units
R = diff(reshape(axis(AX),2,2));    %Gets the axis data range

set(AX,'Units','inches');   %converts axis to absolute units
P = get(AX,'Position');     %Gets current axis position
if R(1) > R(2)  %If the horizontal domain is larger than the vertical, shrink the vertical
%     set(AX,'Position',[P(1:3) P(3)*R(2)/R(1)]);
    W = P(3)*R(2)/R(1);
    set(AX,'Position',[P(1) P(2)+(P(4)-W)/2 P(3) W]);
else    %Else, shrink the horizontal
    W = P(4)*R(1)/R(2);
    set(AX,'Position',[P(1)+(P(3)-W)/2 P(2) W P(4)]);
end

set(AX,'Units',UN); %convert axis back to original units


