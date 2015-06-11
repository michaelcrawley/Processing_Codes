function [varargout] = LES_PhaseAverage(trigger,varargin)

    %Find length of partial forcing period at end of signal
    N = length(trigger);
    Nfp = unique(diff(find(diff(trigger) > 0))); %number of points in forcing period
    cut = mod(N,Nfp); %remove partial period at end of set
    nf = (N-cut)/Nfp;
    
    varargout = cell(size(varargin));
    for n = 1:length(varargin)
        %find which dimension is time
        NV = size(varargin{n});
        ndim = length(NV);
        t_indx = (N == NV)*(1:ndim)';
        indx = repmat({':'},ndim,1);
        indx{t_indx} = 1:NV(t_indx)-cut;
        
        %Reorder for permutation (so that reshape works properly)
        order = [setdiff(1:ndim,t_indx), t_indx];
        
        %Phase average, inverse permute
        tmp = mean(reshape(permute(varargin{n}(indx{:}),order),[NV(order(1:end-1)),Nfp,nf]),ndim+1);
        varargout{n} = ipermute(tmp,order);
    end

end