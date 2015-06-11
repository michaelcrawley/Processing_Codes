function [vm,rms_vel] =horizontal_vrms(z,y,u,lm);

%calculating the horizontal RMS velocity
rms_vel=zeros(1,size(z,2));
vm=0;
mxz=size(z,1);
mnz=1;
for l = 1 : size(z,1)
    if z(l,1)>lm && z(l-1,1)<lm
        mxz=l;
        break
    end
end
for l = 1 : size(z,1)
        if z(l,1)<-lm && z(l+1,1)>-lm
        mnz=l;
        break
    end
end

for j=2:size(z,2)
    count=0;
    for i=mnz:mxz
        rms_vel(j)=rms_vel(j)+(u(i,j))^2;
        count=count+1;
    end
rms_vel(j)=sqrt(rms_vel(j)./count);
    
end

vm=max(rms_vel);
