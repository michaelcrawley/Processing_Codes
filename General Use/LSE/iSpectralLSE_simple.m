function [recon] = iSpectralLSE_simple(c,Ns)

    N = size(c);
    if mod(N(1),2) %odd initial sample length - nyquist frequency not evaluated
        c_up = zeros((N(1)-1)*Ns+1,N(2:end));
        c_up(1:(N(1)-1)/2+1,:,:) = c(1:(N(1)-1)/2+1,:,:);
        c_up(end-(N(1)-1)/2+1:end,:,:) = c((N(1)-1)/2+2:end,:,:);  
    else %even initial sample length
        c_up = zeros(N(1)*Ns,N(2:end));
        c(N(1)/2+1,:,:) = c(N(1)/2+1,:,:)/2; %halve energy at nyquist frequency, since we are placing it twice after upsampling
        c_up(1:N(1)/2+1,:,:) = c(1:N(1)/2+1,:,:);
        c_up(end-(N(1)/2)+1:end,:,:) = c(N(1)/2+1:end,:,:);        
    end    
    c_up = c_up*Ns;
    
    recon = ifft(c_up,[],1);
end