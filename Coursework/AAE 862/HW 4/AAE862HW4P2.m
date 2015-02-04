function [pm pa px] = AAE862HW4P2(U2,yw,alpha1,alpha2)
    %Completed by Michael Crawley
    %Calc mass averaged pressure loss
    Ue = max(U2);
    temp1 = (U2/Ue).*(1-(U2/Ue).^2);
    temp2 = (1-(U2/Ue).^2);
    temp3 = (U2/Ue).*(1-(U2/Ue));
    dl = yw(2:end)-yw(1:end-1);
    pm = (cosd(alpha1)/cosd(alpha2))^2*temp1(2:end)'*dl;
    pa = (cosd(alpha1)/cosd(alpha2))^2*temp2(2:end)'*dl;
    px = 2*temp3(2:end)'*dl;
end