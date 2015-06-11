function [estimates, fCurve, sse] = FSSfit(xdata, ydata)
%This function performs a fit to the Fine Scale Similarity spectrum (see 
%below) using the Nelder-Mead Simplex Method as the core operation. The 
%Nelder-Mead Method is then used to find the best fit x- and y- offset. A 
%similar approach could be used to fit many different functions.

start_point = rand(1, 2);
model = @FSSfun;
estimates = fminsearch(model, start_point);

[sse, fCurve] = model(estimates);

% LSSfun accepts curve offsets as inputs and outputs sse,
% the sum of squares error, and the FittedCurve. FMINSEARCH only needs sse, 
% but we want to plot the FittedCurve at the end.
    function [sse, FittedCurve] = FSSfun(params)
        B = params(1);
        fF = params(2);
        
        fofF = xdata/fF;

        F = zeros(size(ydata));

        F(fofF <= 0.05) = 9.9+14.91126*log10(fofF(fofF <= 0.05));

        F((fofF <= 0.15)&(fofF > 0.05)) = -3.5 + (11.874876 + 2.1202444*log10(20/3*fofF((fofF <= 0.15)&(fofF > 0.05)))+...
            7.5211814*log10(20/3*fofF((fofF <= 0.15)&(fofF > 0.05))).^2).*log10(20/3*fofF((fofF <= 0.15)&(fofF > 0.05)));

        F((fofF <= 1)&(fofF > 0.15)) = (-1.0550362 + 4.9774046*log10(fofF((fofF <= 1)&(fofF > 0.15)))).*log10(fofF((fofF <= 1)&(fofF > 0.15))).^2;

        F((fofF <= 10)&(fofF > 1)) = -(8.1476823 + 3.6523177*log10(fofF((fofF <= 10)&(fofF > 1)))).*log10(fofF((fofF <= 10)&(fofF > 1))).^2;

        F((fofF <= 30)&(fofF > 10)) = -11.8 - (27.2523 + 0.8091863*log10(fofF((fofF <= 30)&(fofF > 10))/10)+...
            14.851964*log10(fofF((fofF <= 30)&(fofF > 10))/10).^2).*log10(fofF((fofF <= 30)&(fofF > 10))/10);

        F(fofF > 30) = 29.77786-38.16739*log10(fofF(fofF > 30));

        FittedCurve = B+F;
        
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector.^2./xdata);
    end
end