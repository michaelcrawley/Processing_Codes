function mn=bl_power_law_fit(y,u,delta,uFS)
%finds the power n that gives a best fit line for a given data set such
%that u/u_freestream=(y/delta)^(1/n)

n=1;
n_min=1;
n_max=20;
mn=1;
total_error_min=inf;
model=1;
error_n=0;
total_error=0;
i=1;

for n=n_min:.01:n_max
    model=uFS.*(y./delta).^(1/n);
    error_n=model-u;
    total_error(i)=sqrt(sum(error_n.^2));
    if total_error(i)<total_error_min
        total_error_min=total_error(i);
        mn=n;
    end
    i=i+1;
end

%plot([n_min:0.01:n_max],total_error);