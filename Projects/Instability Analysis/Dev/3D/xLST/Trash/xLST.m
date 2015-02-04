function [lalpha lphi] = xLST(m,nmode,tr,omega,ialpha,x,y,U)
    %This code performs spatial linear stability theory calculations on PIV
    %data throughout the entire domain.
    [M N] = size(U);
    if (length(x) ~= M && length(x) ~= N) || (length(y) ~= M && length(y) ~= N)
        error('Input matrix size mismatch');
    elseif (length(x) ~= N && length(x) == M)
        U = U';
    end
    
    lalpha(1) = ialpha;
    for i = 1:N
        if i > 2
            lalpha(i) = interp1(x(1:i-1),lalpha(1:i-1),x(i),'spline'); 
        end
        [lalpha(i),lphi{i}] = LSTv2(m,nmode,tr,U(:,i),y,omega,lalpha(end));
        if isnan(lalpha(i))
           warning(['Local Eigenvalue not found at x/D = ',num2str(x(i)),', ending program']);
           break;
        else
            disp(['Processing completed for x/D =',num2str(x(i)),', found eigenvalue = ',num2str(lalpha(i))]);
        end
    end
end