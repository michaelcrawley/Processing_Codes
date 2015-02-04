function [estimates, fCurve, sse] = LSSfit(xdata, ydata)
%This function performs a fit to the Large Scale Similarity spectrum (see 
%below) using the Nelder-Mead Simplex Method as the core operation. The 
%Nelder-Mead Method is then used to find the best fit x- and y- offset. A 
%similar approach could be used to fit many different functions.

start_point = rand(1, 2);
model = @LSSfun;
estimates = fminsearch(model, start_point);

[sse, fCurve] = model(estimates);

% LSSfun accepts curve offsets as inputs and outputs sse,
% the sum of squares error, and the FittedCurve. FMINSEARCH only needs sse, 
% but we want to plot the FittedCurve at the end.
    function [sse, FittedCurve] = LSSfun(params)
        A = params(1);
        fL = params(2);
        
        fofL = xdata/fL;

        L = zeros(size(ydata));

        L(fofL <= 0.5) = 2.53895+18.4*log10(fofL(fofL <= 0.5));

        L((fofL <= 1)&(fofL > 0.5)) = -38.19338*log10(fofL((fofL <= 1)&(fofL > 0.5))).^2-...
            16.91175*log10(fofL((fofL <= 1)&(fofL > 0.5))).^3;

        L((fofL <= 2.5)&(fofL > 1)) = (1.06617 - 45.2994*log10(fofL((fofL <= 2.5)&(fofL > 1)))+...
            21.40972*log10(fofL((fofL <= 2.5)&(fofL > 1))).^2).*log10(fofL((fofL <= 2.5)&(fofL > 1)));

        L(fofL > 2.5) = 5.64174-27.7472*log10(fofL(fofL > 2.5));

        FittedCurve = A+L;
        
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector.^2./xdata);
    end
end