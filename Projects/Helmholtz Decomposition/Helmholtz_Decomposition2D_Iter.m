function [potential,solenoidal,MSE] = Helmholtz_Decomposition2D_Iter(x,y,U,V)

    N = size(U);
    max_itr = 1e4;
    max_error = sqrt(eps);
    
    %Make sure matrices are in correct ascending order
    [~,Ix] = sort(x(:,1));
    [~,Iy] = sort(y(1,:));
    x = x(Ix,:);
    y = y(:,Iy);
    U = U(Ix,Iy);
    V = V(Ix,Iy);
    
    %Create derivative matrices for divergence calculation
    order = 2;
    dx = mean(diff(x(:,1)));
    dy = mean(diff(y(1,:))); %any differences in the spacing will be due to truncation, since the data is assumed to be from PIV
    partial_x_D1 = mNumericalDerivative(1,order,dx,N(1));
    partial_y_D1 = mNumericalDerivative(1,order,dy,N(2));
    
    R = partial_x_D1*V - (partial_y_D1*U')'; %for now we are just going to assume that U and V are 2-D
    
    %initialize output matrices
    solenoidal.u = U;
    solenoidal.v = V;
        
    for n = 1:max_itr
        
        %Divergence updating
        Dt = partial_x_D1*solenoidal.u + (partial_y_D1*solenoidal.v')';
        uD = -0.5*dx*Dt;
        vD = -0.5*dy*Dt;
        
        %Compute updates for interior nodes        
        solenoidal.u(3:end,:) = solenoidal.u(3:end,:) + uD(2:end-1,:);
        solenoidal.u(1:end-2,:) = solenoidal.u(1:end-2,:) - uD(2:end-1,:);
        solenoidal.v(:,3:end) = solenoidal.v(:,3:end) + vD(:,2:end-1);
        solenoidal.v(:,1:end-2) = solenoidal.v(:,1:end-2) - vD(:,2:end-1);
        
        %Curl updating
        Rt = partial_x_D1*solenoidal.v - (partial_y_D1*solenoidal.u')';
        uR = 0.5*dy*(Rt-R);
        vR = -0.5*dx*(Rt-R);
        
        %compute updates for interior nodes
        solenoidal.u(:,3:end) = solenoidal.u(:,3:end) + uR(:,2:end-1);
        solenoidal.u(:,1:end-2) = solenoidal.u(:,1:end-2) - uR(:,2:end-1);
        solenoidal.v(3:end,:) = solenoidal.v(3:end,:) + vR(2:end-1,:);
        solenoidal.v(1:end-2,:) = solenoidal.v(1:end-2,:) - vR(2:end-1,:);
        
        %Compute MSE - mean divergence of solenoidal field
        div = partial_x_D1*solenoidal.u + (partial_y_D1*solenoidal.v')';
        MSE(n) = mean(div(:).^2);
    end
 
    potential.u = U-solenoidal.u;
    potential.v = V-solenoidal.v;
end