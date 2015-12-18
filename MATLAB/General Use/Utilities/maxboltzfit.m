function [params, fCurve, sse] = maxboltzfit(xdata, ydata)
%This function performs a fit to a Maxwell-Boltzmann PDF function (see below) 
%using the Nelder-Mead Simplex Method as the core operation. The range and 
%offset of the data are determined by averaging the beginning and ending 
%10% of ydata. The data is vertically shifted and rescaled to take on the 
%range of -1 to 1. The Nelder-Mead Method is then used to find the best fit
%rate and inflection point. A similar approach could be used to fit many
%different functions.
%
%FUNCTION FORM: y = A*heaviside(x-xo)*(x-xo)^2*exp(-1/2*(x-xo)^2/S^2) +yoff
%   A - Amplitude
%   S - Width
%   xo - centroid
%   yoff - vertical offset
%
%OUTPUT
% Params - contains the four coefficients in this order: .
% fCurve - contains the data points for the fitted curve.
% sse - contains the sum of square errors between ydata and FittedCurve.

L = length(ydata);
RN = mean([mean(ydata(1:ceil(L/100))) mean(ydata(floor(99*L/100):L))]);    %Data Range
yoff = RN; %vertical offset
ydata = ydata - yoff;

[A,I] = max(ydata);
if abs(min(ydata)) > A
    [A,I] = min(ydata);
end
A = mean(ydata(I-ceil(L/100):I+ceil(L/100)));   %Amplitude

ydata = ydata/A;  %scale data for range 0 to 1

    %Finds the remaining curve parameters using 
    %the Nelder-Mead Simplex Method.
start_point = rand(1,3);    %uses random starting point
model = @fun;
est = fminsearch(model, start_point);

[sse, fCurve] = model(est); %Computes best fit curve and error
fCurve = A*fCurve/max(fCurve) +yoff;    %Scales and vertically shifts fitted curve to data range

    %Accumulate fit parameters for output
params(1) = sqrt(est(3)/est(1));  %xo
params(2) = yoff;           %yoff
params(3) = A;              %A
params(4) = 1/sqrt(est(1));       %S

% fun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for the FittedCurve.
    function [sse, FittedCurve] = fun(params)
        C1 = params(1);  
        C2 = params(2);  
        C3 = params(3);

        FittedCurve = heaviside(xdata - sqrt(C3/C1)).* (xdata - sqrt(C3/C1)).*  exp(-1/2*(C1*xdata.^2 -C2*xdata +C3));
        FittedCurve(isnan(FittedCurve)) = 0;
        ErrorVector = FittedCurve/max(FittedCurve) - ydata;
        sse = sum(ErrorVector.^2);
    end
end