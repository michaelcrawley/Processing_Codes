function [ybl,np]=BL_height(z,y,u,vm,percent,vmr,lm)
%limiting z plane range
mxz=size(z,1);
mnz=1;
for l = 1 : size(z,1)
    if z(l,1)>lm & z(l-1,1) < lm
        mxz=l;
        break
    end
end
for l = 1 : size(z,1)
        if z(l,1)<-lm & z(l+1,1) > -lm
        mnz=l;
        break
    end
end
%limiting the y plane
ix = size(y,2);

for k = 1:size(y,2)
    if y(1,k) > 15
        ix= k;
        break
    end
end
%calculate BL Thickness
vm=max(vmr(1:ix));
np=0;
 for i = mnz:mxz
     for j = 2:size(z,2)
         if u(i,j)> percent*vm && u(i,j-1)< percent*vm
             np=np+1;
             ybl(np)=y(i,j-1);
             break
         end
     end
 end