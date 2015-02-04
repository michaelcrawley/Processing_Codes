function [Jacob,R,Q,Delta,lambdaim,lambdaip,lambda2m,lambda2p,lambda2logic] = HW3(data)
    
    %% Build Velocity Gradient Tensor
    Jacob = cell(length(data.y),length(data.x)); S = Jacob; Omega = Jacob; 
    for i = 1:length(data.y)
        velocgrad.dudx(i,:) = NumericalDerivative(1,2,mean(diff(data.x)),data.u(i,:));
        velocgrad.dvdx(i,:) = NumericalDerivative(1,2,mean(diff(data.x)),data.v(i,:));
    end
    for j = 1:length(data.x)
        velocgrad.dudy(:,j) = NumericalDerivative(1,2,mean(diff(data.y)),data.u(:,j)');
        velocgrad.dvdy(:,j) = NumericalDerivative(1,2,mean(diff(data.y)),data.v(:,j)');        
    end

    %% Criterions
    R = zeros(length(data.y),length(data.x)); Q = R; Delta = R; lambda2m = R; lambda2p = R; lambda2logic = R; lambdaim = R; lambdaip = R;
    for i = 1:length(data.y)
        for j = 1:length(data.x)
            Jacob{i,j} = [velocgrad.dudx(i,j) velocgrad.dudy(i,j); velocgrad.dvdx(i,j) velocgrad.dvdy(i,j)];
            R(i,j) = det(Jacob{i,j});
            Q(i,j) = -0.5*(sum(sum(Jacob{i,j}.*Jacob{i,j}')));
            Delta(i,j) = (Q(i,j)/3)^3 + (R(i,j)/2)^2;
            splitter1 = eigs(Jacob{i,j});
            lambdaip(i,j) = splitter1(1,1); lambdaim(i,j) = splitter1(2,1);
            S{i,j}= 0.5*(Jacob{i,j}+Jacob{i,j}');
            Omega{i,j} = 0.5*(Jacob{i,j}-Jacob{i,j}');
            splitter2 = eigs(S{i,j}^2+Omega{i,j}^2);
            lambda2p(i,j) = splitter2(1,1); lambda2m(i,j) = splitter2(2,1);
            if (lambda2m(i,j) < 0)  && (lambda2p(i,j) < 0)
                lambda2logic(i,j) = 1;
            end
        end
    end
    figure;
    subplot(2,1,1);pcolor(data.x,data.y,Delta);title('\Delta Criterion');colorbar;colormap gray; 
    subplot(2,1,2);pcolor(data.x,data.y,Q);title('Q Criterion');colorbar;colormap gray;
    figure;
    subplot(2,1,1);contourf(data.x,data.y,lambdaim);title('\lambda_i_- (Swirling Strength)');colorbar;colormap gray;
    subplot(2,1,2);contourf(data.x,data.y,imag(lambdaim));title('\lambda_i_- Criterion');colorbar; colormap gray;
    figure;
    subplot(2,1,1);contourf(data.x,data.y,lambdaip);title('\lambda_i_+ (Swirling Strength)');colorbar; colormap gray;
    subplot(2,1,2);contourf(data.x,data.y,imag(lambdaip));title('\lambda_i_+ Criterion');colorbar; colormap gray;
    figure;
    subplot(2,1,1);contourf(data.x,data.y,lambda2m);title('\lambda_2_- Criterion');colorbar; colormap gray;
    subplot(2,1,2);contourf(data.x,data.y,lambda2p);title('\lambda_2_+ Criterion');colorbar; colormap gray;
    figure;
    contourf(data.x,data.y,lambda2logic);title('\lambda_2 Criterion (logical)');colorbar; colormap gray; 
end