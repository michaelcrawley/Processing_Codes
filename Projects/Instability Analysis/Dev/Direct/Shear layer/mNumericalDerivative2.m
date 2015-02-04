function [coefm] = mNumericalDerivative2(Dorder, Horder, h, N, varargin)
    %Generate matrix to numerically derivate matrix given constant step size
    %Code Version: 1.0 @ 2011-03-04
    %This is not complete; DRP scheme is only available for central
    %differencing
    %Inputs:
    %   Dorder: Derivative Order
    %   Horder: Number of terms to use
    %   h: step size
    %   N: size of required matrix
    %   options:    '-center' for central difference only
    %               '-upstream' for backward difference only
    %               '-downstream' for forward difference only
    %               '-DRP' for 4th order DRP scheme
    %Outputs:
    %   coefm: matrix for Dorder derivative of phi with Horder accuracy,
    %   Dphi = phi*coefm
    
    %set default scheme (if not given by user)
    if any(strcmp(varargin,'-DRP')) 
        if Dorder ~= 1 || Horder ~= 6
            warning('DRP scheme only available for 1st derivative, 4th order.  Switching to standard Taylor series expansion mode.'); %#ok<WNTAG>
            varargin(strcmp(varargin,'-DRP')) = '';
        end
        if length(varargin) ==1
            varargin{2} = '-center';
        end
    end
    
    if isempty(varargin)
       varargin{1} = '-center'; 
    end
    
    %calculate coefficients for differencing operations
    if any(strcmp(varargin,'-center'))
        if mod(Horder,2)
            error('Given accuracy order cannot be accomplished with central finite difference method; please provide even accuracy order');
        end
        %central differencing coefficients
        n = Horder/2+floor((Dorder-1)/2);
        if any(strcmp(varargin,'-DRP'))
            coefs = DRPlookup(3,3,h);
            coefm = spdiags(repmat(coefs',N,1),-3:3,N,N);
        else
            coefs = TSE(n,n,h,Dorder);
            coefm = spdiags(repmat(coefs',N,1),-n:n,N,N);
        end
        
        %left bound coefficients
        for i = 1:Horder/2+floor((Dorder-1)/2)
            if any(strcmp(varargin,'-DRP'))
                coefm(i,1:7) = DRPlookup(i-1,7-i,h);
            else
                coefm(i,1:2*n+Dorder) = TSE(i-1,Horder+Dorder-i,h,Dorder)';
            end
        end
        
        %right bound coefficients
        for i = N-Horder/2+floor((Dorder-1)/2)+1:N
            if any(strcmp(varargin,'-DRP'))
                coefm(i,end-6:end) = DRPlookup(-N+6+i,N-i,h);
            else
                coefm(i,end+1-(2*n+Dorder):end) = TSE(-N+Horder+Dorder-1+i,N-i,h,Dorder)';
            end
        end
        
    elseif any(strcmp(varargin,'-upstream'));        
        %Upstream differencing coefficients
        coefs = TSE(Horder+Dorder-1,0,h,Dorder);
        coefm = spdiags(repmat(coefs',N,1),-(Horder+Dorder-1):0,N,N);
        
        %Left bound coefficients
        for i = 1:Horder+Dorder-1
           coefm(i,1:Horder+Dorder) = TSE(i-1,Horder+Dorder-i,h,Dorder)';
        end
        
    elseif any(strcmp(varargin,'-downstream'));        
        %Downstream differencing coefficients
        coefs = TSE(0,Horder+Dorder-1,h,Dorder);
        coefm = spdiags(repmat(coefs',N,1),0:Horder+Dorder-1,N,N);
        
        %Right bound coefficients
        for i = N-Horder-Dorder+1:N
           coefm(i,end-(Horder+Dorder-1):end) = TSE(-N+Horder+Dorder-1+i,N-i,h,Dorder)';
        end        
    end
end

function [coefs] = DRPlookup(n,m,h)
%Looks up coefficients for 4th Order, 1st Derivative DRP scheme
    a{33} = [-0.02084314277031176 ...
                0.166705904414580469 ...
                -0.77088238051822552 ...
                0 ...
                0.77088238051822552 ...
                -0.166705904414580469 ...
                0.02084314277031176]';
            
    a{42} = [0.02636943100 ...
                -0.16613853300 ...
                0.518484526d0 ...
                -1.27327473700 ...
                0.47476091400 ...
                0.46884035700 ...
                -0.049041958d0]';
        
    a{51} = [-0.048230454 ...
                0.281814650 ...
                -0.768949766 ...
                1.388928322 ...
                -2.147776050 ...
                1.084875676 ...
                0.209337622]';
            
    a{60} = [0.203876371 ...
                -1.128328861 ...
                2.833498741 ...
                -4.461567104 ...
                5.108851915 ...
                -4.748611401 ...
                2.192280339]';
    
    a{24} = -flipud(a{42});
    a{15} = -flipud(a{51});
    a{06} = -flipud(a{60});
            
    coefs = a{str2double(strcat([num2str(n),num2str(m)]))}/h;
end