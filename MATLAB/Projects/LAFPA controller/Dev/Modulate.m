function [PSt,NPt] = Modulate(F2,Amp,NP,PS,clockpulse)
    NP2 = round(clockpulse/F2);
    NPt = lcm(NP,NP2);
    repnum = NPt/NP; %number of times necessary to rep low frequency pulse sequence
    iAmp = round((Amp/1E6)*clockpulse); %modulation amplitude in indices (Amp variable input is in usec)
    
    [NCh, ~] = size(PS);
    
    %rep the original pulse sequence the required number of times
    PSt = zeros(NCh,2*repnum+2);
    q = PS(:,3) == NP+2; %check if any channels are have a split square wave
    
    %rep non-split waveforms
    for n = 0:repnum-1
        PSt(q,2*n+1:2*n+2) = PS(q,1:2)+n*NP;
        PSt(q,2*n+1:2*n+2) = PSt(q,2*n+1:2*n+2)+repmat(round(iAmp*sin(2*pi*PSt(q,2*n+1)/NP2)),[1 2]);
    end
    PSt(q,end-1:end) = repmat([NPt+1, 1],sum(q),1);
    
    %rep split waveforms
    if ~all(q)
        PS(~q,4) = PS(~q,4)+PS(~q,2)-PS(~q,1)-1;
        PSt(~q,1:2) = PS(~q,1:2);
        for n = 0:repnum-1
            PSt(~q,2*n+3:2*n+4) = PS(~q,3:4)+n*NP;       
            PSt(~q,2*n+3:2*n+4) = PSt(~q,2*n+3:2*n+4)+repmat(round(iAmp*sin(2*pi*PSt(~q,2*n+3)/NP2)),[1 2]);
        end
        q2 = PSt(:,end) > NPt; %test to make sure wavform is still split after modulation
        q3 = ~q & ~q2;
        
        %Modify waveforms that are still split
        PSt(q2,end) = NPt+1;
        PSt(q2,2) = (PSt(q2,4)-PSt(q2,3))-(PSt(q2,end)-PSt(q2,end-1)-1)+1;
        
        %Modify waveforms that are no longer split
        PSt(q3,1:end-2) = PSt(q3,3:end);
        PSt(q3,end-1:end) = repmat([NPt+1 1],sum(q3),1);
    end
       
end