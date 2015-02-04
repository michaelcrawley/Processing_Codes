function [C_delta] = STCWT1D_Cdelta(scales,omegaks,mother,param)
%Calculates the normalization constant C_delta for the (1+1) Dimensional
%Spatio-Temporal Inverse Wavelet Transform. 
%Last Updated: 2014-05-02 by Michael Crawley

    wavelet = MotherWavelets('ST2',mother,param);
    [K,T] = meshgrid(omegaks{2},omegaks{1}); %order is flipped because meshgrid assumes input is DIM2,DIM1
    Wn = zeros(length(scales{1}),length(scales{2}));
    for n = 1:length(scales{1})
        for m = 1:length(scales{2})
            daughter = wavelet(scales{2}(m)*T/sqrt(scales{1}(n)),scales{2}(m)*K*sqrt(scales{1}(n)));
            Wn(n,m) = mean(real(daughter(:)))/scales{1}(n)/scales{2}(m);
        end
    end
    C_delta = trapz(scales{1},trapz(scales{2},Wn,2)); %Integrate over both dimensions
end