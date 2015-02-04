function [PSD f G] = calcPSD2(x,order,dt,whan)

    %Get signal acquisition information
    BS = size(x);    
    
    %Calculate frequency, wavenumber axes
    f = cell(length(order),1);
    for n = 1:length(order)
        nf = ceil((BS(order(n))+1)/2);%number of distinct frequencies returned for PSD
        f{n} = (0:nf-1)'*(1/BS(order(n))/dt(n));   %Frequency axis for dim n
    end
    
    %Create Windows
    wt = whan(BSt); %temporal window
    ws = whan(BSs); %spatial window
    wwt = mean(wt.^2); %window weigth
    wws = mean(ws.^2);
    wt = repmat(wt,[1,BSs,NB]);
    ws = permute(repmat(ws.',[BSt,1,NB]),[2 1 3]);
    
    %Normalize signal
    xm = x - mean(mean(mean(x)));
    
    %Compute FFT
    S = permute(fft(permute(fft(xm.*wt),[2 1 3]).*ws),[2 1 3]);
    
    %Compute one-sided PSD
    G = S(1:nf,1:nk,:); %remove symmetric portion
    PSD = 4*mean(conj(G).*G,3)*(dt/wwt/BSt)*(ds/wws/BSs);
end