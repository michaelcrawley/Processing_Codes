function [Xn Yn varargout] = ReshapeGrid(X,Y,varargin)
%This code sorts matrices based on X and Y locations. Currently, it is
%designed to only be used for the nearfield microphone data, where multiple
%Y locations occur for each X location.

    %Columnize inputs
    X = X(:);
    Y = Y(:);
    for n = 1:(nargin-2)
        if ndims(varargin{n}) > 2
            s = size(varargin{n});
            varargin{n} = reshape(varargin{n},[s(1)*s(2) s(3:end)]);
        else
           varargin{n} = varargin{n}(:); 
        end
    end   
    
    %Initialize new matrices
    uniqueX = unique(X);
    N = length(uniqueX);
    M = length(unique(Y(X == uniqueX(1))));
    Xn = zeros(M,N);
    Yn = Xn;
    varargout = cell(1,nargin-2);
    for n = 1:(nargin-2)
        s = size(varargin{n});
        if iscell(varargin{n})
            varargout{n} = cell(M,N);
        else
            if s(2) ~= 1
                varargout{n} = zeros([M N s(2:end)]);
            else
                varargout{n} = Xn;
            end
        end
    end
    
    %Reorder matrices
    for n = 1:N
        I = X == uniqueX(n);
        [Ys Is] = sort(Y(I));
        
        Xn(:,n) = X(I);
        Yn(:,n) = Ys;
        
       %reorder additional variables
       for nn = 1:(nargin-2)
           CI = repmat({':'},1,ndims(varargin{nn}));
           Cin = CI;
           Cout = repmat({':'},1,ndims(varargout{nn}));
           CI(1) = {I};
           Cin(1) = {Is};
           Cout(2) = {n};

           uniqueV = varargin{nn}(CI{:});
           varargout{nn}(Cout{:}) = uniqueV(Cin{:});
       end
    end
end