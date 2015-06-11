function QQ = PIV_corr(x,v)

%Find correlation Length
q = [];
qq = [];
S = size(v);
for m = 1:S(3)
    progress(m,1,S(3),10);
    for n = 1:S(2)
        q(:,n) = xcorr(v(:,n,m),'coeff');
    end
    qq(:,:,m) = q(S(1):end,:);
end
Q = mean(qq,3);

[C,Ix] = min(Q,[],1);
[C,Iy] = min(C);
Ix = Ix(Iy);

yM = max(find(Q(Ix,:) < mean(Q(Ix,:))))+2;
ym = min(find(Q(Ix,:) < mean(Q(Ix,:))))-2;

QQ = mean(Q(:,ym:yM),2);
% plot(x(:,1),QQ)

