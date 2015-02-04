function [pulselocs A] = CalcPulseLocations(M1,M2,N)

    if ~exist('N','var'), N = 8; end

    IN = (0:N-1);
    SP1 = rem(2*pi*M1*IN/N,2*pi); %mode 1
    SP2 = rem(2*pi*M2*IN/N,2*pi); %mode 2

    SP1(SP1 < 0) = SP1(SP1 < 0) + 2*pi;
    SP2(SP2 < 0) = SP2(SP2 < 0) + 2*pi;
    
    %Phase location of combined waveform
    Q = abs(SP1-SP2) > pi;
    pulselocs = (SP1+SP2)/2; %pulse location in terms of phase

    %if combined signal amplitude < 0.94, then no pulse
    A = zeros(N,1); %amplitude (used to turn off channels)
    for n = 1:N
        if Q(n) == 1
            pulselocs(n) = pulselocs(n) - pi;
            A(n) = 2*pi-abs(SP1(n)-SP2(n)) <= 11*pi/16;
        else
            A(n) = abs(SP1(n)-SP2(n)) <= 11*pi/16;
        end
    end
end