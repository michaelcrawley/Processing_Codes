function [yn,I,BL,DT,MT,H] = BLcalc(x,y)
%For calculating the properties of a shear layer based on a profile
%given by x and y. It calculates the various properties in two ways:
%integrating from 10% to 98%, and linear fitting the steepest point and
%extrapolating that line back to where u=0-then integrating over that
%fitted profile.
%
%INPUTS
% x - abcissa vector
% y - ordinate vector
%
%OUTPUTS
% yn - y is interpolated and curve fit. yn contains the resulting profile
% I - index numbers corresponding to 10% and 98% limits
% BL - width of shear layer based on 2% and 98% limits
% DT - displacement thickness calculated two ways
% MT - momentum thickness calculated two ways
% H - shape factor calculated two ways
%
%Assumes x is decreasing and y is increasing as a function of index number

    %Interpolates data for better determination of threshold points
xn = linspace(x(end),x(1),length(x)*3);
xn = fliplr(xn)';
y = interp1(x,y,xn);
x = xn;

Q = round(length(y)/10);

I = zeros(2,1); %Locates 10% and 98% points
I(1) = min(find(y >= 0.10+0.90*mean(y(1:Q))));
I(2) = min(find(y >= 0.98));

IL = min(find(y >= 0.02+0.98*mean(y(1:Q))));
BL = x(IL)-x(I(2)); %Calculate BL width

dy = diff(y);
[Q,Im] = max(diff(y));  %Find steepest point of profile

C = polyfit(x(Im:Im+2),y(Im:Im+2),1);   %Fit a few points nearest the steepest point
yn = y;

x0 = -C(2)/C(1);
x1 = min(find(x <= x0));    %locates x-intercept of linear fit
yn(x1:Im-1) = C(1)*x(x1:Im-1)+C(2); %create fitted profile
yn(1:x1-1) = 0;

DT = zeros(2,1); MT = DT; H = DT;
DT(1) = -trapz(x(x1-1:I(2)),1-yn(x1-1:I(2))); %Linear fit method
MT(1) = -trapz(x(x1-1:I(2)),yn(x1-1:I(2)).*(1-yn(x1-1:I(2))));

DT(2) = -trapz(x(I(1):I(2)),1-y(I(1):I(2)));   %Simple 10% to 98% method
MT(2) = -trapz(x(I(1):I(2)),y(I(1):I(2)).*(1-y(I(1):I(2))));

H = DT./MT;
