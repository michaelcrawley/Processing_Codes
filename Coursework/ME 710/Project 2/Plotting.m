aloc = [0.001 0.004 0.01 0.04 0.08 0.10 0.20];
figure;hold on;
for i = 1:length(aloc)
    [~, loc] = min(abs(xs-aloc(i)));
    plot(rs,(T(:,loc,2)-Tw)/(Tm(1,loc)-Tw));
end
title('Nondimensional Temperature profile at selected x^* locations');xlabel('r^*');ylabel('\theta');xlim([0 1]);ylim([0 2]);legend(num2str(aloc'));