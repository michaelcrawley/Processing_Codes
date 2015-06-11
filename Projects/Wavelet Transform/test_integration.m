

wavelet = MotherWavelets('ST2','morlet',[6 6]);

%Frequency indices
omegaks{1} = linspace(0,6.2832e5/2,1025);
omegaks{2} = linspace(0,125.0912,33);
% [K,T] = meshgrid(omegaks{2},omegaks{1}); %order is flipped because meshgrid assumes input is DIM2,DIM1

%Scale indices
scales{1} = linspace(1.2263,1.5474e5,100);
scales{2} = linspace(6.8605e-4,0.2484,50);
% [S2,S1] = meshgrid(scales{2},scales{1}); %order is flipped because meshgrid assumes input is DIM2,DIM1

% %Integrate over frequencies first, then scales
% Wn = zeros(length(scales{1}),length(scales{2}));
% for n = 1:length(scales{1})
%     for m = 1:length(scales{2})
%         daughter = wavelet(scales{2}(m)*T/sqrt(scales{1}(n)),scales{2}(m)*K*sqrt(scales{1}(n)));
%         Wn(n,m) = mean(real(daughter(:)))/scales{1}(n)/scales{2}(m);
%     end
% end
% C_delta1 = trapz(scales{1},trapz(scales{2},Wn,2)); %Integrate over both dimensions
% 
% %Integrate over scales first, the frequencies
% Wn = zeros(length(omegaks{1}),length(omegaks{2}));
% for n = 1:length(omegaks{1})
%     for m = 1:length(omegaks{2})
%         daughter = wavelet(S2*omegaks{1}(n)./sqrt(S1),S2*omegaks{2}(m).*sqrt(S1))./S1./S2;
%         Wn(n,m) = mean(real(daughter(:)));
%     end
% end
% C_delta2 = trapz(omegaks{1},trapz(omegaks{2},Wn,2)); %Integrate over both dimensions

wavelets = zeros(length(omegaks{1}),length(omegaks{2}),length(scales{1}),length(scales{2}));
for n = 1:length(omegaks{1})
    for m = 1:length(omegaks{2})
        for k = 1:length(scales{1})
            for p = 1:length(scales{2})
                wavelets(n,m,k,p) = wavelet(scales{2}(p)*omegaks{1}(n)/sqrt(scales{1}(k)),scales{2}(p)*omegaks{2}(m)*sqrt(scales{1}(k)))/scales{1}(k)/scales{2}(p);
            end
        end
    end
end