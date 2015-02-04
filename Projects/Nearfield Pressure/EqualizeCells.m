function varargout = EqualizeCells(varargin)

    for n = 1:length(varargin)
        L = min(cellfun(@length,varargin{n}(:))); %minimum length
        shortened = cellfun(@(x) x((length(x)-L)/2+1:end-(length(x)-L)/2), varargin{n},'uniformoutput',false);
        dim_cell = (size(shortened{1}) > 1).*(1:length(size(shortened{1})));
        dim_cell = dim_cell(dim_cell > 0);
        ndim_array = max((size(shortened) > 1).*(1:length(size(shortened))));
        dims = 1:ndim_array+length(dim_cell);
        rotation = [setdiff(dims,dim_cell) intersect(dims,dim_cell)];
        rotated = cellfun(@(x) permute(x,rotation),shortened,'uniformoutput',false);
        varargout{n} = cell2mat(rotated);
    end
end