function [Psi PsiT] = orthowavelet_basis(BS)
    %Inputs: BS - block size (number of samples)
    %Outputs:   Psi - LMB wavelet basis (time domain)
    %           PsiT - LMB wavelet basis (Fourier domain)
    %
    %       This program automatically uses LMB wavelet with n = 4
    %       Based on Mallat 1989 algorithm, appendix A
    
   
    
    %Cubic Splines
    N1 = @(omega) 5 + 30*cos(omega/2).^2+30*(sin(omega/2).^2).*(cos(omega/2).^2); 
    N2 = @(omega) 2*(sin(omega/2).^4).*(cos(omega/2).^2) + 70*(cos(omega/2).^4) + 2/3*(sin(omega/2).^6);
    Sigma8 = @(omega) (N1(omega)+N2(omega))./(105*(sin(omega/2).^8));
    
    omega = linspace(-20*pi,20*pi,BS); %frequency resolutions
    
    PsiT = exp(-1i*omega/2).*sqrt(Sigma8(omega/2+pi))./(omega.^4)./sqrt(Sigma8(omega).*Sigma8(omega/2));
    PsiT = ifftshift(PsiT);
    omega = ifftshift(omega);
    PsiT(isnan(PsiT)) = 0; %fix singularity at n*pi

    Psi = real(ifftshift(ifft(PsiT))); %imaginary part is only due to numerical error
end