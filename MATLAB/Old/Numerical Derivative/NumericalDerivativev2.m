function [deriv] = NumericalDerivativev2(Dorder, Horder, h, phi,varargin)
    %Numerically derivate matrix given constant step size
    %Code Version: 2.0 @ 2010-11-14
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
    deriv = NaN*ones(M,N);
    
    coefs = UnDetCoefs(Horder,h,'center',Dorder);
    coefsL = UnDetCoefs(Horder+Dorder-1,h,'down',Dorder);
    coefsR = UnDetCoefs(Horder+Dorder-1,h,'up',Dorder);
    
    if ~isempty(varargin)
        if strcmp(varargin,'-center')
            for j = 1:M
               for i = Horder+1:N-Horder
                   deriv(j,i) = phi(j,i-Horder:i+Horder)*coefs;
               end
            end
        elseif strcmp(varargin,'-downstream')
            for j = 1:M
               for i = 1:N-Horder-Dorder+1
                   deriv(j,i) = phi(j,i:i+Horder+Dorder-1)*coefsL;
               end
            end            
        elseif strcmp(varargin,'-upstream')
            for j = 1:M
               for i = Horder+Dorder:N
                   deriv(j,i) = phi(j,i-Horder-Dorder+1:i)*coefsR;
               end
            end
        end
    else
        for j = 1:M
            for i = Horder+1:N-Horder
               deriv(j,i) = phi(j,i-Horder:i+Horder)*coefs;
            end
            for i = 1:Horder%Left Bounds
               deriv(j,i) = phi(j,i:i+Horder+Dorder-1)*coefsL; 
            end
            for i = N-Horder+1:N%Right Bounds
               deriv(j,i) = phi(j,i-Horder-Dorder+1:i)*coefsR; 
            end
        end
    end
end