%Reynolds Stress sources - assume coordinates are x/D and y/D

uu = U(:,:,M(M>0)).^2;
uv = U(:,:,M(M>0)).*V(:,:,M(M>0));
vv = V(:,:,M(M>0)).^2;

n = 1;

tmp = zeros(size(U(:,:,n))); dtmpdx = tmp; dtmpdy = tmp;
SRC = zeros(size(uu));

for n = 1:size(uu,3)
    %uu11
    tmp = uu(:,:,n);
    dtmpdx(2:end-1,:,:) = (tmp(3:end,:,:)-tmp(1:end-2,:,:))./(x(3:end,:,:)-x(1:end-2,:,:));
    tmp = dtmpdx;
    dtmpdx(2:end-1,:,:) = (tmp(3:end,:,:)-tmp(1:end-2,:,:))./(x(3:end,:,:)-x(1:end-2,:,:));
    SRC(:,:,n) = SRC(:,:,n) +dtmpdx*(1/Uj)^2;

    %uv12
    tmp = uv(:,:,n);
    dtmpdx(2:end-1,:,:) = (tmp(3:end,:,:)-tmp(1:end-2,:,:))./(x(3:end,:,:)-x(1:end-2,:,:));
    tmp = dtmpdx;
    dtmpdy(:,2:end-1,:) = (tmp(:,3:end,:)-tmp(:,1:end-2,:))./(y(:,3:end,:)-y(:,1:end-2,:));
    SRC(:,:,n) = SRC(:,:,n) +dtmpdy*2*(1/Uj)^2;

    %vv22
    tmp = vv(:,:,n);
    dtmpdy(:,2:end-1,:) = (tmp(:,3:end,:)-tmp(:,1:end-2,:))./(y(:,3:end,:)-y(:,1:end-2,:));
    tmp = dtmpdy;
    dtmpdy(:,2:end-1,:) = (tmp(:,3:end,:)-tmp(:,1:end-2,:))./(y(:,3:end,:)-y(:,1:end-2,:));
    SRC(:,:,n) = SRC(:,:,n) +dtmpdy*(1/Uj)^2;
end
clear tmp dtmpdx dtmpdy n


dtmpdx(2:end-1,:,:) = (tmp(3:end,:,:)-tmp(1:end-2,:,:))./(x(3:end,:,:)-x(1:end-2,:,:));
dtmpdy(:,2:end-1,:) = (tmp(:,3:end,:)-tmp(:,1:end-2,:))./(y(:,3:end,:)-y(:,1:end-2,:));


