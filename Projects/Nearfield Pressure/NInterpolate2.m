function [output] = NInterpolate2(varargin)

    if ~mod(nargin,2), error('Error: Wrong Number of Inputs'); end
    
    N = (nargin-1)/2+1;
    input = varargin{N};
    icoords = varargin(1:N-1); %input data coordinates
    ocoords = varargin(N+1:end); %desired output coordinates
    clear('varargin'); %for memory usage
    
    Lo = size(ocoords{1});
    Li = size(icoords{1});
    S = length(icoords);
    
    output = zeros(numel(ocoords{1}),1);
    for n = 1:numel(ocoords{1})
        
        %Calc distance from desired point to all points
        d = zeros(Li);
        for q = 1:S
            d = d + (ocoords{q}(n)-icoords{q}).^2;
        end
        d = sqrt(d);
        
        %Find index of nearest point
        [~,I] = min(d(:));
        sub = cell(1,length(Li));
        [sub{:}] = ind2sub(Li,I);
        
        %Find direction of next nearest point along each dimension
        indbox = cell(1,S);
        dx = zeros(1,S);
        for q = 1:S
            dx(q) = ocoords{q}(n)-icoords{q}(sub{:}); %distance from desired to nearest, along dim q
            indbox{q} = [sub{q} sub{q}+sign(dx(q))];
        end

        %Find all values in N-D box around desired location and
        %cooresponding indices
        box_vals = cell(1,S+1);
        box_vals{1} = input(indbox{:});
        box_locs = cell(1,S+1);
        for q = 1:S
            box_locs{1}{q} = icoords{q}(indbox{:});
        end 
        
        %Interpolate along all dimensions
        str = repmat({':'},1,S);
        for q = 1:S
            dim = str;
            dim{q} = 1;
            
            %Find new box vertices
            L = diff(box_locs{q}{q},[],q);
            for qq = q:S
                slope = diff(box_locs{q}{qq},[],q)./L;
                slope(isnan(slope)) = 0;
                box_locs{q+1}{qq} = box_locs{q}{qq}(dim{:})+dx(q)*slope;
            end
            
            %Calc slopes and new box vertice values
            slope = diff(box_vals{q},[],q)./L;
            slope(isnan(slope)) = 0;             
            box_vals{q+1} = box_vals{q}(dim{:})+dx(q)*slope;
        end
        
        output(n) = squeeze(box_vals{end});
    end
    
    output = reshape(output,Lo);
    
end