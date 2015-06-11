function [] = Project705()
%Completed by Michael Crawley 2009-12-04 as part of the requirements for ME
%705

%Input values
CarLength = 180*0.0254; %in meters
CarHeight = 52*0.0254; %in meters
CarWidth = 60*0.0254; %in meters
CarVelocity = 75*0.44704; %in meters per second
xlocation = (0.5*CarLength-3*0.0254); %x location of spoiler offset (z plane) in meters
ylocation = 33*0.0254; %y location of spoiler offset (z plane) in meters
l = 6.25*.0254; %cord length of airfoil in meters
t = 2*.0254; %maximum thickness of airfoil in meters
alpha=-7;
density = 1.177; %air density at standard atmospheric conditions

%calculate ellipse flow field
x1 = 0:0.1:5;
y1 = 0:0.1:2;
a1 = 0.5*(CarHeight+0.5*CarLength);
c1 = (a1*(0.5*CarLength-a1))^0.5;

u1=zeros(length(x1),length(y1));
v1=zeros(length(x1),length(y1));

for i = 1:length(x1)
    for j = 1:length(y1)
        [u1(i,j) v1(i,j)]=EllipseFlowField(CarVelocity,a1,c1,x1(i),y1(j));
    end    
end

xe=0:.0001:0.5*CarLength;
ellipse=zeros(1,length(xe));
for i=1:length(xe)
    ellipse(i)=sqrt(1-(xe(i)/(0.5*CarLength))^2)*CarHeight;
end

figure
hold on
quiver(x1,y1,u1',v1')
plot(xlocation,ylocation,'o r')
plot(xe,ellipse,'g')
hold off

%calculate velocity components at airfoil location in ellipse flow field
[ue ve]=EllipseFlowField(CarVelocity,a1,c1,xlocation,ylocation);
R = sqrt(ue^2+ve^2);
theta = (180/pi)*atan(ve/ue);

%calculate lift force
UnitForce = density*pi*(R^2)*(l+0.77*t)*sin((theta+alpha)*pi/180);
TotalForce = UnitForce*CarWidth*0.224808943; %in pounds
fprintf('The force developed is (lbs): %6.2f\n',TotalForce);

%Parametric Studies
x=xlocation-6*.0254:0.01:xlocation+12*0.0254;
y=ylocation-2*.0254:0.01:ylocation+12*0.0254;

for i=1:length(x)
    [ux(i) vx(i)]=EllipseFlowField(CarVelocity,a1,c1,x(i),ylocation);
    Rx(i) = sqrt(ux(i)^2+vx(i)^2);
    thetax(i) = (180/pi)*atan(vx(i)/ux(i));
    UnitForcex(i) = density*pi*(Rx(i)^2)*(l+0.77*t)*sin((thetax(i)+alpha)*pi/180);
    TotalForcex(i) = UnitForcex(i)*CarWidth*0.224808943; 
end

for j=1:length(y)
    [uy(j) vy(j)]=EllipseFlowField(CarVelocity,a1,c1,xlocation,y(j));
    Ry(j) = sqrt(uy(j)^2+vy(j)^2);
    thetay(j) = (180/pi)*atan(vy(j)/uy(j));
    UnitForcey(j) = density*pi*(Ry(j)^2)*(l+0.77*t)*sin((thetay(j)+alpha)*pi/180);
    TotalForcey(j) = UnitForcey(j)*CarWidth*0.224808943; 
end

figure
plot(x-xlocation,-TotalForcex);title('Downforce as a function of length (m)');xlabel('Distance Behind Original Placement (m)');ylabel('Downforce (lbs)');
figure
plot(y-ylocation,-TotalForcey);title('Downforce as a function of height (m)');xlabel('Distance Above Original Placement (m)');ylabel('Downforce (lbs)');
end

function [u,v]=EllipseFlowField(U,a,c,x,y)
    z=complex(x,y);
    W=U*(1+((a/c)^2-1)*(0.5-z/(4*sqrt((z/2)^2-c^2))));
    u=real(W);
    v=-imag(W);
end