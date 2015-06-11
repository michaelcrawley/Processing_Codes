function [xn] = iSTCWT1D(wave,scales,omegaks,mother,param,C_delta)
%Computes (1+1)-Dimensional Inverse Spatio-Temporal Continuous Wavelet
%Transform using Fourier method. Based off of Torrence 1998. 
%Inputs:
%           wave:       wavelet (Q,J,M,N), where DIM3 is time
%           scales:     scale vectors (cell array)
%           omegaks:    Fourier dimensions (cell array, optional)
%           mother:     mother wavelet (Default is Morlet) (optional)
%           param:      parameter for mother wavelet (2,1) (optional)
%           C_delta:    Normalization energy of the mother wavelet
%                       (optional)
%Outputs:
%           xn:         Reconstructed signal (M,N)
%Last Updated: 2014-05-02 by Michael Crawley

    %Integrate over scales - take real part only
    wave = real(wave);
    N = size(wave);
    [a,c] = meshgrid(scales{2},scales{1}); %meshgrid assumes inputs are DIM2,DIM1
    a = repmat(a,[1 1 N(3:4)]);
    c = repmat(c,[1 1 N(3:4)]);
    wave = wave./a./a./c; %normalize by spatial scale and phase velocity
    xn = trapz(scales{2},wave,2); %integrate along spatial scales
    xn = trapz(scales{1},xn,1); %integrate along phase velocities
    xn = squeeze(xn);
    
    %Calculate C_delta if not provided
    if ~exist('C_delta','var')||isempty(C_delta)
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
    
    %Normalize reconstructed signal
    xn = xn/C_delta;
end