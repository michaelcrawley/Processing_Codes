function [estparams, FittedCurve, newfun,sse] = GeneralFit(x,y,fun,xo)
%Inputs: x data, y data, Function handle, initial guess (optional)
%Format for function needs to be '@(x,etc)' where x is the domain variable,
%and variable names cannot be the letter 'o'

    %Rotate arrays if necessary
    if  (sum((size(x) == size(y))) == 0) && (sum((size(x) == size(y'))) == 2)
        y = y';
    elseif (sum((size(x) == size(y))) <= 1) && (sum((size(x) == size(y'))) <= 1)
        error('Data array size mismatch');
    end

    %Find parameters
    funstr = func2str(fun);
    lastchar = find(funstr == ')',1);
    parameters = length(find(funstr(1:lastchar) == ','));
    locations = [2 find((funstr(1:lastchar) == ',')) lastchar];
    newfun = funstr(lastchar+1:end);
    for i = 2:length(locations)-1
       newfun =  strrep(newfun,funstr(locations(i)+1:locations(i+1)-1),['o(',num2str(i-1),')']);
    end
    newfun = strcat('@(x,o)',newfun);
    newfun = str2func(newfun);    
    
    %Call FMINSEARCH using random start location for parameters
    if ~exist('xo','var'), xo = rand(1,parameters); end
    estparams = fminsearch(@generalfun,xo);   
    [FittedCurve] = newfun(x,estparams);
    sse = sum((FittedCurve-y).^2);
    
    function [sse, CurveFit] = generalfun(params)
        CurveFit = newfun(x,params);
        sse = sum((CurveFit - y).^2);        
    end
end