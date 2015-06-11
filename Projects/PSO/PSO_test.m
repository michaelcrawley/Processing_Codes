tic;
fun = @(x) 10*x.^3-2*x.^2+5*x+1;
N = 1e7;
imax = 1e2;

chi = 0.5;
alpha = 1;
c1 = 0.5;
c2 = 0.5;
phi1 = rand(N,imax);
phi2 = rand(N,imax);

position = 100*(rand(N,1)-0.5);
velocity = 50*(rand(N,1)-0.5);

value = fun(position);
individual_best.position = position;
individual_best.value = value;
[global_best.value,I] = min(abs(value));
global_best.position = position(I);

for n = 1:imax
    velocity = chi*(alpha*velocity + c1*phi1(:,n).*(individual_best.position-position) + c2*phi2(:,n).*(global_best.position-position));
    position = position + velocity;
    
    value = fun(position);
    
    individual_best.position = position.*(abs(value) <= individual_best.value) + individual_best.position.*(abs(value) > individual_best.value);
    individual_best.value = value.*(abs(value) <= individual_best.value) + individual_best.value.*(abs(value) > individual_best.value);
    
    [current_best.value,I] = min(abs(value));
    current_best.position = position(I);
    [global_best.value,I] = min(abs([global_best.value,current_best.value]));
    global_best.position = global_best.position*(I==1) + current_best.position*(I==2);    
end
toc