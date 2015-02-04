function [] = Project705Previous()
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
theta = (180/pi)*atan(ue/ve);

%calculate airfoil flow field
x2=0:.01:2*l;
y2=-2*t:.01:2*t;
a2=0.25*(l+0.77*t);
c2=l/4;

u2=zeros(length(x2),length(y2));
v2=zeros(length(x2),length(y2));

for i=1:length(x2)
    for j=1:length(y2)
        [u2(i,j) v2(i,j)]=AirfoilFlowField(R,a2,c2,t,theta+alpha,x2(i),y2(j));
    end
end

figure
quiver(x2,y2,u2',v2')

%calculate velocity components in the ellipse flow field versus the airfoil
%flow field
[ua va]=EllipseFlowField(CarVelocity,a1,c1,xlocation,ylocation - 3*0.0254);
[ub vb]=AirfoilFlowField(CarVelocity,a2,c2,t,theta+alpha,0,-3*0.0254);
fprintf('Velocity components in the Ellipse flow field: u=%f, v=%f\n',ua,va);
fprintf('Velocity components in the Airfoil flow field: u=%f, v=%f\n',ub,vb);

%calculate lift force
UnitForce = density*pi*(R^2)*(l+0.77*t)*sin((theta+alpha)*pi/180);
TotalForce = UnitForce*CarWidth*0.224808943; %in pounds
fprintf('The force developed is (lbs): %6.2f\n',TotalForce);

%Parametric Studies
x=xlocation:0.01:xlocation+12*0.0254;
y=ylocation:0.01:ylocation+12*0.0254;

for i=1:length(x)
    [ux(i) vx(i)]=EllipseFlowField(CarVelocity,a1,c1,x(i),ylocation);
    Rx(i) = sqrt(ux(i)^2+vx(i)^2);
    thetax(i) = (180/pi)*atan(ux(i)/vx(i));
    UnitForcex(i) = density*pi*(Rx(i)^2)*(l+0.77*t)*sin((thetax(i)+alpha)*pi/180);
    TotalForcex(i) = UnitForcex(i)*CarWidth*0.224808943; 
end

for j=1:length(y)
    [uy(j) vy(j)]=EllipseFlowField(CarVelocity,a1,c1,xlocation,y(j));
    Ry(j) = sqrt(uy(j)^2+vy(j)^2);
    thetay(j) = (180/pi)*atan(uy(j)/vy(j));
    UnitForcey(j) = density*pi*(Ry(j)^2)*(l+0.77*t)*sin((thetay(j)+alpha)*pi/180);
    TotalForcey(j) = UnitForcey(j)*CarWidth*0.224808943; 
end

figure
plot(x-xlocation,-TotalForcex);title('Downforce as a function of length(m)');
figure
plot(y-ylocation,-TotalForcey);title('Downforce as a function of Height(m)');

end

function [u,v]=EllipseFlowField(U,a,c,x,y)
    z=complex(x,y);
    W=U*(1+((a/c)^2-1)*(0.5-z/(4*sqrt((z/2)^2-c^2))));
    u=real(W);
    v=-imag(W);
end

function [u,v]=AirfoilFlowField(U,a,c,t,theta,x,y)
    z=complex(x,y);
    gamma=4*pi*U*a*sin(theta*pi/180);
    temp1=0.5+z/(4*sqrt((z/2)^2-c^2));
    temp2=z/2+sqrt((z/2)^2-c^2)+(0.77*t/4);
    W=U*(temp1*exp(complex(0,-theta))-(a^2)*temp1/(temp2^2)*exp(complex(0,theta)))+complex(0,gamma/(2*pi))*(temp1/temp2);
    u=real(W);
    v=-imag(W);
end