function [PW] = calcPW(Freq,Duty) 
    dfreq = Freq/1000;
    if Duty < 0
        if dfreq>30
            tf1 = 0.2857;
            tf2 = 11.4286;
        else
            tf1 = 0.6;
            tf2 = 2;
        end
        PW = ((((dfreq*tf1)+tf2)*10)/dfreq)/1000000;
    else
        per = 1/Freq;
        PW = Duty*per/100;
    end    
end