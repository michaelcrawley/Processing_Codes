function [PA,PAm] = calcPA(PS, N, NP,CH,RForDC)
    PAm = zeros(N,NP);
    converter = 2.^((1:N)-1);
    
    
    for i = 1:N
        counter = 1;
        qswitch = RForDC;
        if CH(i) == 0
            PAm(i,:) = RForDC;
        else 
            for j = 1:NP
                if j == PS(i,counter)
                    qswitch = ~qswitch;
                    counter = counter+1;
                end
                PAm(i,j) = qswitch;
            end  
        end
    end
    PA = converter*PAm;
end