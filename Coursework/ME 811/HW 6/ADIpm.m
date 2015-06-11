function [phi R] = ADIpm(A,Rtol,itrmax)
    [M N] = size(A.O);
    
    if ~exist('Rtol','var')
        Rtol = 1E-5;
    end
    if ~exist('itrmax','var')
        itrmax = 1E4;
    end 
    
    phi = zeros(M,N);
    R = zeros(1,itrmax);
    [~,R(1)] = calcRes(A,phi);
    counter = 0;
    while (R(counter+1) >= Rtol) && (counter < itrmax)
        %Row sweep
        S = A.P(1,:) - A.N(1,:).*phi(2,:);
		phi(1,:) = TDMAsolver(A.W(1,:),A.O(1,:),A.E(1,:),S);
		for m=2:M-1
			S = A.P(m,:) - A.S(m,:).*phi(m-1,:) - A.N(m,:).*phi(m+1,:);
			phi(m,:) = TDMAsolver(A.W(m,:),A.O(m,:),A.E(m,:),S);
		end
		S = A.P(M,:) - phi(M-1,:).*A.S(M,:);
		phi(M,:) = TDMAsolver(A.W(M,:),A.O(M,:),A.E(M,:),S);
        
        %Column sweep
        S = A.P(:,1) - A.E(:,1).*phi(:,2);
		phi(:,1) = TDMAsolver(A.S(:,1),A.O(:,1),A.N(:,1),S);
		for n=2:N-1
			S = A.P(:,n) - A.W(:,n).*phi(:,n-1) - A.E(:,n).*phi(:,n+1);
			phi(:,n) = TDMAsolver(A.S(:,n),A.O(:,n),A.N(:,n),S);
		end
		S = A.P(:,N) - A.W(:,N).*phi(:,N-1);
		phi(:,N) = TDMAsolver(A.S(:,N),A.O(:,N),A.N(:,N),S);        
        
        counter = counter + 1;
        [~,R(counter+1)] = calcRes(A,phi);
    end
    R = R(1:counter+1);
end