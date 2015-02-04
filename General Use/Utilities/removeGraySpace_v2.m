function removeGraySpace_v2(FIG,varargin)
%
%Program removes useless gray space around axes including axes labels. 
%Allows user to specify the amount of space they want between rows and 
%columns. Program assumes that figure was produced by subplot commands with
%no subsequent tinkering. Also, program cannot handle subplots with merged
%elements. To avoid statistical discrepancies for certain combinations,
%position origin coordinates are rounded to the nearest 5%. If you require
%more precision, change the parameter PR immediately below this help
%section, but know that certain combinations of axes will not be cropped 
%properly. All desired figure elements, except a global title for subplots, 
%should be in place before executing this program. Global titles should be
%added after execution by the appropriate manipulation of figure units and
%sizes. Annotated figures with dissassociated text boxes are clearly beyond
%the scope of this program and should be added after gray space removal.
%However, the command "text" produces objects which will move with the
%appropriate axis.
%
%GENERAL RULE: Do not manipulate figure with mouse before running this
%program - Undesirable results may occur due to broken relationships
%between figure objects.
%
%Program accepts the following input types:
%   removeGraySpace(FIG) - This sets spacing between rows/columns to zero.
%       This is also the proper call for a figure with only one axis.
%   removeGraySpace(FIG,rSpace,cSpace) - Uses specified spacing values.
%
%FIG - handle to desired figure (e.g. gcf)
%cSpace (percentage) - Spacing between columns beyond necessary minimum
%rSpace (percentage) - Spacing between rows beyond necessary minimum
%

LE = 0.01;  %Left/Right Edge Padding (percentage)
TE = 0.01;  %Top/Bottom Edge Padding (percentage)
PR = 0.05;  %Position precision (percentage)

if ~isempty(varargin)
    rSpace = varargin{1};
    cSpace = varargin{2};
else
    cSpace = 0;
    rSpace = 0;
end

Q = get(FIG,'Children');
if isempty(Q)
    error('No Axis handles found!');    %Spits an error if you try to pass a bad figure handle
end
A = [];
C = [];
for n = length(Q):-1:1
    if isempty(strfind(get(Q(n),'Type'),'ui'))
        T = get(Q(n),'Tag');
        if strcmp(T,'Colorbar') %look for objects which are colorbars
            C(end+1) = Q(n);
        elseif strcmp(T,'legend')   %look for objects which are legends
            C(end+1) = Q(n);
        else
            A(end+1) = Q(n);
        end
    else
        Q = [Q(1:n-1) Q(n+1:end)];
    end 
end

FU = get(FIG,'Units');  %Original figure units
AU = [];    %Original axes units cell array
for n = 1:length(A)
    AU{n} = get(A(n),'Units');
    set(A(n),'Units','normalized'); %convert to normalized units
end

OP = [];
for n = 1:length(A)
    OP(n,:) = get(A(n),'OuterPosition');
end
X = [OP(:,1:2) OP(:,1:2)+OP(:,3:4)];
AX = zeros(size(C));
CP = [];
LOC = [];
for n = 1:length(C) %Find appropriate axis for non-axis objects
    CP(n,:) = get(C(n),'Position');
    CTI(n,:) = get(C(n),'TightInset');
    LOC{n} = get(C(n),'Location');
    Q = (CP(n,1) < X(:,3)) & (CP(n,1) > X(:,1));
    Q = Q & (CP(n,2) < X(:,4)) & (CP(n,2) > X(:,2));
    AX(n) = find(Q);
end

%Position => [in from left, in from bottom, width, height]
for n = 1:length(A)
    P(n,:) = get(A(n),'Position');  %Stores position of axes
    TI(n,:) = get(A(n),'TightInset');   %Stores needed margin around axes for labels, etc.
end
P(:,1:2) = round(P(:,1:2)/PR)*PR;   %Prevents small discrepancies from derailing program

for n = 1:length(C)
    Q = findstr('Outside',LOC{n});  %If the object is outside the graph
    if ~isempty(Q)
        Q = LOC{n}(1:Q-1);
        if strcmp(Q,'North')    %Add space to TightInset on appropriate side for object        
            TI(AX(n),4) = CP(n,4) +CTI(n,4)+0.03 +CP(n,2) -(P(AX(n),2)+P(AX(n),4));
        elseif strcmp(Q,'East')
            TI(AX(n),3) = TI(AX(n),3) +CP(n,3) +CTI(n,3)+0.03 +CP(n,1) -(P(AX(n),1)+P(AX(n),3));
        elseif strcmp(Q,'South')
            TI(AX(n),2) = CTI(n,2) +P(AX(n),2) -CP(n,2);
        else
            TI(AX(n),1) = CTI(n,1) +P(AX(n),1) -CP(n,1);
        end
    end
end

    %Moving axes require changing the units of the figure. Original
    %units will be restored before program terminates. 

set(FIG,'Units','inches');    
set(A(1),'Units','inches');
Q = get(A(1),'Position');
set(A(1),'Units','normalized');
DU = Q(3)/P(1,3);   %conversion for percentage to inches
DU(2) = Q(4)/P(1,4);


%%%% ADJUST COLUMNS %%%%    
S0 = 0; %Axis column offset from previous column
S1 = 0; %Column width
LM = 0;  %Axis column offset from left figure margin
UN = unique(P(:,1));    %Find locations of left margins of columns

for m = 1:length(UN)    %Iterate through each column
    X = find(P(:,1) == UN(m));
    for n = 1:length(X) %Iterate through rows in mth column looking for maximum needed offset.
        if TI(X(n),1) > S0
            S0 = TI(X(n),1);
        end
    end
    for n = 1:length(X) %Look for largest needed width for mth column
        if P(X(n),3)+S0+TI(X(n),3) > S1
            S1 = P(X(n),3)+S0+TI(X(n),3);
        end
    end
    if m==1
        P(X,1) = LM +S0 +LE;  %Overwrite left offset of first row in mth column
        LM = LM +LE;
    else
        P(X,1) = LM +S0;  %Overwrite left offset of other rows in mth column
    end
    LM = LM +S1 +cSpace +0.01;  %calculate LM offset for (m+1)th column
    S0 = 0;  %reset maximum column offset
    S1 = 0;
end


%%%% ADJUST ROWS %%%%
S0 = 0; %Axis row offset from previous row
S1 = 0; %row height
BM = 0;  %Axis height offset from bottom figure margin
UN = unique(P(:,2));    %Find locations of bottom margins of rows

for m = 1:length(UN)    %Iterate through each row
    X = find(P(:,2) == UN(m));
    for n = 1:length(X) %Iterate through columns in mth column looking for maximum needed offset.
        if TI(X(n),2) > S0
            S0 = TI(X(n),2);
        end
    end
    for n = 1:length(X) %Look for largest needed height for mth row
        if P(X(n),4)+S0+TI(X(n),4) > S1
            S1 = P(X(n),4)+S0+TI(X(n),4);
        end
    end
    if m==1
        P(X,2) = BM +S0 +TE;  %Overwrite bottom offset of first column in mth row
        BM = BM +TE;
    else
        P(X,2) = BM +S0;  %Overwrite bottom offset of other columns in mth row
    end
    BM = BM +S1 +rSpace +0.01;  %calculate BM offset for (m+1)th row
    S0 = 0;  %reset maximum row offset
    S1 = 0;
end

for n = 1:length(A) %Apply new position coordinates
    set(A(n),'Position',P(n,:));
    set(A(n),'Units','inches'); %convert axes to inches so that they won't scale when figure is shrunk to wrap them
end

OP = get(FIG,'OuterPosition');  %Resize figure to wrap tightly around total needed area
set(FIG,'OuterPosition',[OP(1) OP(2) (LM-cSpace+2*LE)*DU(1) 0.8+(BM-rSpace+2*TE)*DU(2)]);
for n = 1:length(A)
    set(A(n),'Units','normalized'); %convert to normalized units
end

    %Restore original units for all objects
for n = 1:length(A)
    set(A(n),'Units',AU{n});
end
set(FIG,'Units',FU);