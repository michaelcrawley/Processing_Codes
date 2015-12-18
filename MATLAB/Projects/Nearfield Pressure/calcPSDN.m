function [PSD f S] = calcPSDN(x,order,dt,whan)

    %Get processing information
    N = length(order);

    %Get signal acquisition information
    BS = size(x);  
    NBS = length(BS);
    
    %Calculate frequency, wavenumber axes
    f = cell(N,1);
    for n = 1:N
        nf = ceil((BS(order(n))+1)/2);%number of distinct frequencies returned for PSD
        f{n} = (0:nf-1)'*(1/BS(order(n))/dt(n));   %Frequency axis for dim n
    end
    
    %Create Windows
    wndo = cell(N,1); %full window matrix
    wt = zeros(N,1); %window weight
    for n = 1:N
        tmp = whan(BS(order(n))); %single window
        wt(n) = mean(tmp.^2);
        tmp = permute(tmp,circshift(1:NBS,[1 order(n)-1]));
        reps = BS; 
        reps(order(n)) = 1;
        wndo{n} = repmat(tmp,reps);
    end
        
    %Normalize signal
    xm = x;
    for n = 1:N
        xm = mean(xm,order(n));        
    end
    u = setdiff(1:NBS,order);
    reps = BS;
    for n = 1:length(u)
        reps(u(n)) = 1;        
    end
    xm = repmat(xm,reps);
    xn = x-xm;
    
    %Compute FFT
    S = xn;
    for n = 1:N
        S = fft(S.*wndo{n},[],order(n));
    end
    
    %Compute single-sided PSD
    PSD = (2^N)*(conj(S).*S)*prod(dt)/prod(wt)/prod(BS);
end