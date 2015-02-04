function [deriv] = NumericalDerivative(Dorder, Horder, h, phi)
    %Numerically derivate array given constant step size
    %Inputs:
    %Dorder: Derivative Order
    %Horder: Accuracy Order
    %h: step size
    %phi: function to derivate
    %Outputs:
    %deriv: Dorder derivative of phi with Horder accuracy
    [M N] = size(phi);
    
    coefs = UnDetCoefs(-Horder,Horder,h,Dorder);
    coefsL = UnDetCoefs(0,Horder+1,h,Dorder);
    coefsR = UnDetCoefs(-Horder-1,0,h,Dorder);
    
    deriv = zeros(size(phi));
    for j = 1:M
        for i = Horder+1:N-Horder
           deriv(j,i) = phi(j,i-Horder:i+Horder)*coefs;
        end
        for i = 1:Horder%Left Bounds
           deriv(j,i) = phi(j,i:i+Horder+1)*coefsL; 
        end
        for i = N-Horder+1:N%Right Bounds
           deriv(j,i) = phi(j,i-Horder-1:i)*coefsR; 
        end
    end
end