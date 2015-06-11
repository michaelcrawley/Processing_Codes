%Problem 2, Homework 5, ME 710

Lx = [ 0:10:190 0:10:190 zeros(1,5) 200*ones(1,5) ] ;
Rx = [10:10:200 10:10:200 zeros(1,5) 200*ones(1,5)];
Ty = [zeros(1,20) 50*ones(1,20) 10:10:50 10:10:50 ];
Ly = [zeros(1,20) 50*ones(1,20) 0:10:40 0:10:40 ];

F = zeros(length(Lx));
for i = 1:length(Lx)
    for j = 1:length(Lx)
        if i ~= j
            d1 = sqrt((Rx(j)-Lx(i))^2+(Ly(i)-Ty(j))^2);
            d2 = sqrt((Lx(j)-Rx(i))^2+(Ly(j)-Ty(i))^2);
            s1 = sqrt((Rx(j)-Rx(i))^2+(Ty(j)-Ty(i))^2);
            s2 = sqrt((Lx(j)-Lx(i))^2+(Ly(j)-Ly(i))^2);
            F(i,j) = abs(d1+d2-s1-s2)/20;
        end
    end
end

sum(F)

%Heating element is located at F(11)
%Slab is located at F(21:40)

plot(F(11,21:40));