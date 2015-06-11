function [deriv] = NumericalDerivativev3(Dorder, Horder, h, phi,varargin)
    %Numerically derivate matrix given constant step size
    %Code Version: 3.0 @ 2011-02-14
    %Inputs:
    %   Dorder: Derivative Order
    %   Horder: Number of terms to use
    %   h: step size
    %   phi: function to derivate
    %   options:    '-center' for central difference only
    %               '-upstream' for backward difference only
    %               '-downstream' for forward difference only
    %Outputs:
    %   deriv: Dorder derivative of phi with Horder accuracy
    
    [M N] = size(phi);
    deriv = zeros(M,N);
    
    %set default scheme (if not given by user)
    if isempty(varargin)
       varargin{1} = '-center'; 
    end
    
    %perform differencing operations
    if strcmp(varargin,'-center')
        if mod(Horder,2)
            error('Given accuracy order cannot be accomplished with central finite difference method; please provide even accuracy order');
        end
        %central differencing
        coefs = TSE(Horder/2+floor((Dorder-1)/2),Horder/2+floor((Dorder-1)/2),h,Dorder);
        for i = 1+Horder/2+floor((Dorder-1)/2):N-Horder/2+floor((Dorder-1)/2)
            deriv(:,i) = phi(:,i-Horder/2+floor((Dorder-1)/2):i+Horder/2+floor((Dorder-1)/2))*coefs;
        end
        
        %left bound 
        for i = 1:Horder/2+floor((Dorder-1)/2)
           coefs = TSE(i-1,Horder+Dorder-i,h,Dorder);
           deriv(:,i) = phi(:,1:Horder+Dorder)*coefs; 
        end
        
        %right bound 
        for i = N-Horder/2+floor((Dorder-1)/2)+1:N
           coefs = TSE(-N+Horder+Dorder-1+i,N-i,h,Dorder);
           deriv(:,i) = phi(:,N-Horder-Dorder+1:end)*coefs; 
        end
    elseif strcmp(varargin,'-upstream');        
        %Upstream differencing
        coefs = TSE(Horder+Dorder-1,0,h,Dorder);
        for i = Horder+Dorder:N
          deriv(:,i) = phi(:,i-Horder-Dorder+1:i)*coefs; 
        end
        
        %Left bound
        for i = 1:Horder+Dorder-1
           coefs = TSE(i-1,Horder+Dorder-i,h,Dorder);
           deriv(:,i) = phi(:,1:Horder+Dorder)*coefs; 
        end
    elseif strcmp(varargin,'-downstream');        
        %Downstream differencing
        coefs = TSE(0,Horder+Dorder-1,h,Dorder);
        for i = 1:N-Horder-Dorder+1
          deriv(:,i) = phi(:,i:i+Horder+Dorder-1)*coefs; 
        end  
        
        %Right bound
        for i = N-Horder-Dorder:N
           coefs = TSE(-N+Horder+Dorder-1+i,N-i,h,Dorder);
           deriv(:,i) = phi(:,N-Horder-Dorder+1:end)*coefs; 
        end
    end
end