function [output] = NInterpolate(varargin)

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
    parfor n = 1:numel(ocoords{1})
        
        %Calc distance from desired point to all points
        d = zeros(Li);
        for q = 1:S
            d = d + (ocoords{q}(n)-icoords{q}).^2;
        end
        d = sqrt(d);
        
        %Find index of nearest point
        [~,I] = min(d(:));
        sub = cell(1,length(Li));
        for q = S:-1:2
            k = floor((I-1)/prod(Li(1:q-1)));
            sub{q} = k+1;
            I = I-k*prod(Li(1:q-1));
        end
        sub{1} = I;
        nearest_val = input(sub{:});
        
        %Calc change based off of slope at nearest point
        chn = 0;
        for q = 1:S
            dx = ocoords{q}(n)-icoords{q}(sub{:}); %distance from desired to nearest, along dim q
            
            if dx %only perform calculations if desired point is off of dimension (otherwise delx = NaN)
                subm = sub;
                subp = sub;
                if sub{q} == Li(q)
                    loc_dir = false;                    
                    subm{q} = subm{q}-1; 
                elseif sub{q} == 1
                    loc_dir = true;
                    subp{q} = subp{q}+1;                    
                else
                    subm{q} = subm{q}-1;
                    subp{q} = subp{q}+1;
                    loc_dir = d(subp{q}) < d(subm{q});
                end                

                if loc_dir %desired is closer in positive direction
                    nsub = subp;
                else %desired point is closer in negative direction
                    nsub = subm;
                end

                delx = (input(nsub{:})-input(sub{:}))/(icoords{q}(nsub{:})-icoords{q}(sub{:}));
                chn = chn + delx*dx;
            end
        end
        output(n) = nearest_val + chn;
        
    end
    
    output = reshape(output,Lo);
    
end