D = 0.0254;
Uj = 288;
x = z(1,41:end-10)/D;
y = 1.20 + (x-x(1))*tand(8.6);
T = 2048;

L = length(x);
nearfield = zeros(T,L);
counter = 1;
for n = 1:L
    [~,indx] = min(abs(y(n)-r(:,1)/D));
    nearfield(:,counter) = squeeze(ps(indx,41+n,1:T));
    counter = counter +1;
end


