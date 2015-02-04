function [] = H2P1()
    dxyr = 1/1023;
    yr = 0:dxyr:1;       
    Me = sqrt(1-0.75^2);
    EQNH2P1 = @(yr) (902.56.*(yr.^5) - 875.52.*(yr.^4) + 290.72.*(yr.^3) - 45.727.*(yr.^2) + 4.5034.*(yr) + 0.1454).*((0 <= yr) & (yr < 0.4)) + (1733.3.*(yr.^6) - 5782.1.*(yr.^5) + 7955.3.*(yr.^4) - 5787.3.*(yr.^3) + 2350.5.*(yr.^2) - 507.09.*(yr) + 45.666).*((0.4 <= yr) & (yr <= 0.85)) + (-21.*(yr.^2) + 42.35.*(yr) - 21.198).*((0.85 < yr) & (yr <= 1.0));
    Decay = @(k,xr) exp(2*pi*(k-1)*xr);
    Decay2 = @(k,xr) exp(2*pi*(k-1)*xr/Me);
    Pprime = EQNH2P1(yr);    
    Amp = max(Pprime) - min(Pprime);
    FFT_vals = fft(Pprime);
    xr = [0 -.05 -.1 -.2 -.5];
    k = 1:length(yr);    
    
    error = 1;
    n = -1;
    while(error > .01 && n<length(yr))
        n = n+1;
        temp = FFT_vals;
        temp(n+1:end-n) = 0;
        rebuild = ifft(temp);
        error = sqrt(sum((Pprime-rebuild).^2)/length(Pprime))/Amp;
    end
    
    DecayExp = zeros(length(xr),length(yr));
    Upstream = zeros(size(DecayExp));
    for i = 1:length(xr)
        DecayExp(i,:) = Decay(k,xr(i));
        Upstream(i,:) = real(ifft(DecayExp(i,:).*temp));        
    end
    
    error = 1;
    xri = 0;
    dxr = .001;
    while error > .01 
        xri = xri - dxr;
        DecayExpi = Decay(k,xri);
        Upstreami = real(ifft(DecayExpi.*temp));
        error = (max(Upstreami)-min(Upstreami))/Amp;
    end

    error = 1;
    xrb = 0;
    dxrb = .001;
    while error > .01 
        xrb = xrb - dxrb;
        DecayExpb = Decay2(k,xrb);
        Upstreamb = real(ifft(DecayExpb.*temp));
        error = (max(Upstreamb)-min(Upstreamb))/Amp;
    end    
    
               
    fprintf('To account for 99 percent of the energy, only %d wavenumbers are needed\n', n);
    fprintf('Minimum upstream distance where pressure perturbation is 1 percent: %1.3f\n', xri);
    fprintf('Minimum upstream distance where pressure perturbation is 1 percent with a background flow of Mach 0.75: %1.3f\n', xrb);
    figure; plot(yr, Pprime, yr, real(rebuild)); xlabel('y/w'); ylabel('P'); title('True and Reduced Wavenumber Pressure Fields'); legend('True Pressure Field', 'Reduced Pressure Field');
    figure; plot(abs(FFT_vals)); xlabel('wavenumber'); ylabel('amplitude'); xlim([1 1001]);
    figure; plot(abs(FFT_vals)); xlabel('wavenumber'); ylabel('amplitude'); xlim([1 20]);
    figure; plot(yr, Upstream(1,:), yr, Upstream(2,:), yr, Upstream(3,:), yr, Upstream(4,:), yr, Upstream(5,:));  xlabel('y/w'); ylabel('P'); title('Upstream Pressure Fields'); legend('0', '-w/20', '-w/10', '-w/5', '-w/2');
end