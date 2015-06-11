function AAE862HW4P1(PS,SS,cd,alpha,Rec)
    %Completed By Michael Crawley for AAE 862 HW#4 Problem 1
    %Part A
    SurfLen = @(x,y) sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
    PS.dL = SurfLen(PS.x,PS.y);
    SS.dL = SurfLen(SS.x,SS.y);
    
    PS.L = zeros(length(PS.cp),1);
    SS.L = zeros(length(SS.cp),1);
    PS.pl = zeros(length(PS.cp),1);
    SS.pl = zeros(length(SS.cp),1);
    
    W = PS.x(end)/0.99;
    for i = 1:length(PS.L)-1
        PS.L(i+1) = PS.L(i)+PS.dL(i);
        PS.pl(i+1) = PS.pl(i)+(PS.cp(i)^(3/2))*PS.dL(i);
    end    
    for i = 1:length(SS.L)-1
        SS.L(i+1) = SS.L(i)+SS.dL(i);
        SS.pl(i+1) = SS.pl(i)+(SS.cp(i)^(3/2))*SS.dL(i);
    end
    
    PS.plT = PS.pl*2*cd/(W*cos(alpha*pi/180));
    SS.plT = SS.pl*2*cd/(W*cos(alpha*pi/180));
    
    figure;plot(PS.x,PS.y,'k',SS.x,SS.y,'k','LineWidth',2);title('Blade profile');
    figure;plot(PS.L,PS.cp,'k',SS.L,SS.cp,'--k','LineWidth',2);title('Cp Dist vs surface length (Turbulent BL)');xlabel('Distance (mm)');ylabel('Pressure Coef');legend('Pressure Surface','Suction Surface','Location','Best');
    figure;plot(PS.L,PS.plT,'k',SS.L,SS.plT,'--k','LineWidth',2);title('Total Pressure Loss vs. surface length (Turbulent BL)'); xlabel('Distance (mm)');ylabel('Normalized Pressure loss');legend('Pressure Surface','Suction Surface','Location','Best');

    %Part C
    PS.cpx = zeros(length(PS.cp),1);
    SS.cpx = zeros(length(SS.cp),1);
    PS.cd = zeros(length(PS.cp),1);
    SS.cd = zeros(length(SS.cp),1);
    SS.plLi = zeros(length(SS.cp), length(Rec));
    PS.plLi = zeros(length(PS.cp), length(Rec));
    
    for i=1:length(PS.L)-1
        PS.cpx(i+1) = PS.cpx(i) + (PS.cp(i)^(5/2))*PS.dL(i);
    end
    for i=1:length(SS.L)-1
        SS.cpx(i+1) = SS.cpx(i) + (SS.cp(i)^(5/2))*SS.dL(i);
    end
    
    for j = 1:length(Rec)
        for i = 1:length(SS.L)
            if SS.cpx(i) ~= 0
                SS.Ret(i,j) = sqrt(0.45*Rec(j)*SS.cpx(i)/(SS.x(end)*SS.cp(i)^2));
                SS.cd(i,j) = 0.173/SS.Ret(i,j); 
            else
                SS.Ret(i,j) = 0;
                SS.cd(i,j) = 0;
            end
        end
        for i = 1:length(PS.L)
            if PS.cpx(i) ~= 0 
                PS.Ret(i,j) = sqrt(0.45*Rec(j)*PS.cpx(i)/(PS.x(end)*PS.cp(i)^2));
                PS.cd(i,j) = 0.173/PS.Ret(i,j);
            else
                PS.Ret(i,j) = 0;
                PS.cd(i,j) = 0;
            end
        end
        
        for i = 1:length(SS.L) - 1
            SS.plLi(i+1,j) = SS.plLi(i,j)+(SS.cp(i)^(3/2)*SS.cd(i,j))*SS.dL(i);
        end
        for i = 1:length(PS.L) - 1
            PS.plLi(i+1,j) = PS.plLi(i,j)+(PS.cp(i)^(3/2)*PS.cd(i,j))*PS.dL(i);
        end
        SS.plL(:,j)= SS.plLi(:,j)*2/(W*cos(alpha*pi/180));
        PS.plL(:,j)= PS.plLi(:,j)*2/(W*cos(alpha*pi/180));
        figure;plot(PS.L,PS.Ret(:,j),'k',SS.L,SS.Ret(:,j),'--k','LineWidth',2);title({['Re\theta vs. surface length (Laminar BL) Reynolds Number: ' num2str(Rec(j),'%.0f')]});xlabel('Distance (mm)');ylabel('Re\theta');legend('Pressure Surface','Suction Surface','Location','Best');
        figure;plot(PS.L,PS.cd(:,j),'k',SS.L,SS.cd(:,j),'--k','LineWidth',2);title({['C_d vs. surface length (Laminar BL) Reynolds Number: ' num2str(Rec(j),'%.0f')]});xlabel('Distance (mm)');ylabel('C_d');legend('Pressure Surface','Suction Surface','Location','Best');
        figure;plot(PS.L,PS.plL(:,j),'k',SS.L,SS.plL(:,j),'--k','LineWidth',2);title({['Total Pressure Loss vs. surface length (Laminar BL) Reynolds Number: ' num2str(Rec(j),'%.0f')]}); xlabel('Distance (mm)');ylabel('Normalized Pressure loss');legend('Pressure Surface','Suction Surface','Location','Best');
    end
end