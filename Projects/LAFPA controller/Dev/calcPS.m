function [PST, NP] = calcPS(F,D,M,M2,N,rotFlag,RA,CH,clockpulse)

    NP = round(clockpulse/F);
    PW = calcPW(F,D);
    %N = 8;      %Number of channels
    %M = 1;      %Azimuthal mode
    %M2 = -1;   %Second azimuthal mode

          %starting phase for each channel
    IN = zeros(N,1);
    IN(1:N) = (0:N-1);
    SP = rem(2*pi*M*IN/N,2*pi);    %mode 1
    SP2 = rem(2*pi*M2*IN/N,2*pi);     %mode 2

    AN = zeros(N,1);
    for n = 1:N
        AN(n) = n;
        if SP(n) < 0
             SP(n) = SP(n) + 2*pi;
        end
        if SP2(n) < 0
             SP2(n) = SP2(n) + 2*pi;
        end
    end

         %phase location of combined waveform
    Q = abs(SP-SP2) > pi;
    P = (SP+SP2)/2;

         %if combined signal amplitude < 0.94, then no pulse
    A = zeros(N,1);
    for n = 1:N
        if Q(n)==1
            P(n) = P(n) - pi;
            A(n) = 2*pi - abs(SP(n)-SP2(n)) <= 11*pi/16;
        else
            A(n) = abs(SP(n)-SP2(n)) <= 11*pi/16;
        end
    end

        %Find state switching locations - NP - Pulses per period,
        %F - forcing frequency, PW - pulse width

    PS = zeros(N,4);
    PS(1:N,1) = round(P/2/pi*NP);   %index number of begin of pulse
    PS(1:N,2) = PS(1:N,1) + round(PW*F*NP);   


    if rotFlag==1 
    C=NR*N;  %calculate total number of pulse sequences in rotational period
    X=2*C+2;
    PST = zeros(N, X);
    n=zeros(N,1);
    periodcounter=0;

    for i=1:N
        for j=1:NR
            for k=1:N
                S=k+(i-1);
                if S>N
                    S=rem(S,N);
                end
                if A(S)~=0
                    PST(k,1+2*n(k))=PS(S,1)+periodcounter*NP;
                    PST(k,2+2*n(k))=PS(S,2)+periodcounter*NP;
                    n(k)=n(k)+1;
                end
            end
            periodcounter=periodcounter+1;
        end
    end

    for l=1:N
        PST(l,1+2*n(k))=C*NP+1; %set end value
        if PS(l,2) > NP %Account for split wave forms
            for m=1:X-2
                PST(l+1,X+1-m)=PST(l+1,X-m-1); %Shift data
            end
            PST(l+1,1)=0; %initialize first cell
            PST(l+1,2)=PS(l,2)-NP; %initialize second cell
        end
    end
    PS = PST;


    else
        C=1;
        X = 4;
        for n = 1:N
            if PS(n,2) > NP
                PS(n,1:4) = [0 PS(n,2)-NP PS(n,1) NP+1];
            end
            if PS(n,3)==0
                PS(n,3) = NP+1;
            end
            if A(n)==0
                PS(n,1) = NP+1;
            end
        end
    PST=PS;
    end

    AN = circshift(AN,RA);
    for p=1:N
        PST(p,1:X)=PS(AN(p),1:X);
        if CH(p)==0
            PST(p, 1)=C*NP+1;
        end
    end
    
    PST = PST+1;
end
