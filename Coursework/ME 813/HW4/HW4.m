function [B,phi,lambda,lambda_test,modes,a,ac,Error,eigenfunction_check] = HW4(data)
    %Code completed by Michael Crawley for ME 813 HW#4

    B = (data.u*data.u')/length(data.t);figure;contourf(B);title('Cross-correlation tensor');xlabel('x');ylabel('x');
    scale = 1/(length(data.x)-1);
    [phi eigenvalue] = svd(scale*B);
    lambda = diag(eigenvalue);
    EnormT = sum(scale*sum(data.u.*data.u,2))/length(data.t);
    modes = find(cumsum(lambda) >= 0.99*EnormT,1);
    
    figure;hold on;title('Retained Eigenfunction Modes');xlabel('x');
    for i = 1:modes
        plot(data.x,phi(:,i),'Color',rand(3,1));
    end
    hold off;
    
    figure;hold on;title('Modal Coefficients');xlabel('t');
    a = zeros(modes,length(data.t));
    lambda_test = zeros(modes,1);
    for i = 1:modes
        a(i,:) = phi(:,i)'*data.u;
        lambda_test(i,1) = mean(0.04*a(i,:).*a(i,:));
        plot(data.t,a(i,:),'Color',rand(3,1));
    end
    hold off;
    eigenfunction_check = phi'*phi;
    eigenfunction_check(eigenfunction_check < 10*eps) = 0;
    
    I = find(abs(data.t - 3.43) < 10*eps);
    figure;
    subplot(3,1,1);plot(data.x,data.u(:,I));title('True Temperature Distribution');xlabel('x');
    subplot(3,1,2);plot(data.x,phi(:,1:modes)*a(:,I));title('Reconstructed Temperature Distribution');xlabel('x');
    subplot(3,1,3);plot(data.x,data.u(:,I)-phi(:,1:modes)*a(:,I));title('Error in Reconstruction');xlabel('x');
    
    X = [mean(data.u(1,2:end).^2) mean(data.u(1,2:end).*data.u(1,1:end-1)); mean(data.u(1,2:end).*data.u(1,1:end-1)) mean(data.u(1,1:end-1).^2)];
    ac = zeros(modes,length(data.t)-1);
    for i = 1:modes
         Z = [mean(data.u(1,2:end).*a(i,2:end)); mean(data.u(1,1:end-1).*a(i,2:end))];
         D = X\Z;
         ac(i,:) = D(1)*data.u(1,2:end)+D(2)*data.u(1,1:end-1);
    end
    figure;
    subplot(2,1,1);plot(data.t,a(1,:),data.t(2:end),ac(1,:));title('a(1)');xlabel('t');legend('POD Coefficient','SE Coefficient');
    subplot(2,1,2);plot(data.t,a(2,:),data.t(2:end),ac(2,:));title('a(2)');xlabel('t');legend('POD Coefficient','SE Coefficient');
   
    Error = mean((ac-a(:,2:end)).^2,2);
end