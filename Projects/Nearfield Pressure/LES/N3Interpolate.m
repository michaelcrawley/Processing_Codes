function [output] = N3Interpolate(varargin)
%2-D or 3-D interpolation based off of subgrid. Uses TriScatteredInterp
%algorithm. For 2-D interp, inputs are X,Y,V,intX,intY; for 3-D,
%X,Y,Z,V,intX,intY,intZ. In all cases, the dimensions of X,Y,Z,V need to
%match.

    if nargin ~= 5 && nargin ~= 7, 
        error('Error: Wrong Number of Inputs'); 
    end    
    
    L = 1; %length of tesselation grid (+/-L)
    N = (nargin-1)/2+1;
    input = varargin{N};
    icoords = varargin(1:N-1); %input data coordinates
    ocoords = varargin(N+1:end); %desired output coordinates
    clear('varargin'); %for memory usage
    
    Lo = size(ocoords{1});
    Li = size(icoords{1});
    S = length(icoords);
    
    output = zeros(numel(ocoords{1}),1);
    parfor n = 1:numel(ocoords{1})
        
        %Calc distance from desired point to all points
        d = zeros(Li);
        loc = cell(1,3);
        for q = 1:S
            d = d + (ocoords{q}(n)-icoords{q}).^2;
            loc{q} = ocoords{q}(n);
        end
        d = sqrt(d);
        
        %Find index of nearest point
        [~,I] = min(d(:));
        sub = cell(1,length(Li));
        [sub{:}] = ind2sub(Li,I);
        
        %Find SubGrid
        subgridi = cellfun(@(x) (x-L:x+L),sub,'uniformoutput',false); %subgrid indices
        for q = 1:S %get rid of out of bounds grid points
            subgridi{q} = subgridi{q}(subgridi{q} >= 1 & subgridi{q} <= Li(q));
        end
        subgrid = cell(1,S); %subgrid loc vals
        for q = 1:S
            subgrid{q} = icoords{q}(subgridi{:});
            subgrid{q} = subgrid{q}(:);
        end    
        
        %Create tesselation, and interpolate
        vals = input(subgridi{:});
        vals = vals(:);
        F = TriScatteredInterp(subgrid{:},vals);
        output(n) = F(loc{:});
    end
    
    output = reshape(output,Lo);
    
end