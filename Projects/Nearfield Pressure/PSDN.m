function [PSD f S xm wndo] = PSDN(x,order,dt,whan,flag)

    if ~exist('flag','var'), flag = true; end
    
    %Get processing information
    N = length(order);

    %Get signal acquisition information
    BS = size(x);  
    NBS = length(BS);
    
    %Calculate frequency, wavenumber axes
    f = cell(N,1);
    for n = 1:N
        nf = ceil((BS(order(n))+1)/2);%number of distinct frequencies returned for PSD
        if rem(BS(order(n)),2) %NFFT is odd
            f{n} = ifftshift((-(nf-1):nf-1)'*(1/BS(order(n))/dt(n)));   %Full frequency axis for dim n
        else  %NFFT is even          
            f{n} = ifftshift((-(nf-1):nf-2)'*(1/BS(order(n))/dt(n)));   %Full frequency axis for dim n
        end
    end
    
    %Create Window
    wndo = ones(size(x)); %full window matrix
    wt = zeros(N,1); %window weight
    for n = 1:N
        wnd1d = whan(BS(order(n))); %single window
        wt(n) = mean(wnd1d.^2);
        wnd1d = permute(wnd1d,circshift(1:NBS,[1 order(n)-1]));
        reps = BS; 
        reps(order(n)) = 1;
        wnd1d = repmat(wnd1d,reps);
        wndo = wndo.*wnd1d;
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
    if flag %subtract out mean if flag is true
        xn = x-xm;
    else
        xn = x;
    end
    
    %Compute N-Dim FFT
    S = xn.*wndo; %apply window
    for n = 1:N
        S = fft(S,[],order(n));
    end
    
    %Compute single-sided PSD
    PSD = (conj(S).*S)*prod(dt)/mean(wt)/prod(BS(order));
end