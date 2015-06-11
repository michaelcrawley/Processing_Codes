function [ds,corr,c,vn,vp] = schlierenCorrelate(X,Y,W,phD)
% For correlating a segment of a schlieren image (moving window
% correlation). Used in conjunction with schlierenVelocity.m
c = 0;

w = [200, round(W*1.2)];    %pixels - window dimensions (height, width)
s = 0.1;     %i/pixels - window roll-off rate
S = size(X);    %pixels - image size
[y,x] = meshgrid((1:S(2)),(1:S(1)));

xo = [476, 400];   %pixel - window center (vertical, horizontal)
wndw = (tanh(s*(x-xo(1)+w(1)/2))/2 +tanh(-s*(x-xo(1)-w(1)/2))/2) .*(tanh(s*(y-xo(2)+w(2)/2))/2 +tanh(-s*(y-xo(2)-w(2)/2))/2);

I = wndw(:) > 0.01;
IL = min(x(I)); IL(2) = max(x(I)); clear I;
X = X(IL(1):IL(2),:); Y = Y(IL(1):IL(2),:); S = size(X);
[y,x] = meshgrid((1:S(2)),(1:S(1)));

sp = round(w(2)/2);
cs = round(w(2)/20);    %window step size
corr = zeros(ceil((S(2)-2*sp+1)/cs),2*S(2)-1); 
ds = corr; j = 0;
for n = sp:cs:S(2)-sp
    j = j+1;
    progress(n,sp,S(2)-sp,10);
    xo = [S(1)/2, n];   %pixel - window center (vertical, horizontal)
    wndw = (tanh(s*(x-xo(1)+w(1)/2))/2 +tanh(-s*(x-xo(1)-w(1)/2))/2) .*(tanh(s*(y-xo(2)+w(2)/2))/2 +tanh(-s*(y-xo(2)-w(2)/2))/2);
    corr(j,:) = xcorr2_1d(Y,X.*wndw,2);
    ds(j,:) = n-(1-S(2):S(2)-1);
    
    if j==60
        dumb = 0;
    end
end

%%%% PEAK LOCATE
vn = zeros(j,2); vp = vn;
for n = 1:j
    w = moving_average(corr(n,:),3);
    
    vn(n,1) = ds(n,S(2)); vp(n,1) = vn(n,1);
    
    q = round(S(2)+W*phD*[0.65 1.35]);
    [M,I] = max(w(q(1):q(2))); 
    if (I==1) || (I==diff(q)+1)
        vp(n,2) = 0;
    else
        vp(n,2) = q(1)-S(2)+I-1;
    end
    
    q = round(S(2)-W*(1-phD)*[1.35 0.65]);
    [M,I] = max(w(q(1):q(2))); 
    if (I==1) || (I==diff(q)+1)
        vn(n,2) = 0;
    else
        vn(n,2) = q(1)-S(2)+I-1;
    end
    
% % %     q = [0 sign(diff(w))];
% % %     q = ([0 sign(diff(q))] < 0).*(corr(n,:) > 0.2);
% % %     
% % %     in = q==1;
% % %     c(n,1:length(find(in))) = sort(ds(n,in)-ds(n,1280));
% % %     
% % %     q = find(sort(ds(n,in)-ds(n,1280)) < 0,1,'last');
% % %     if ~isempty(q)
% % %         vn(n,:) = [ds(n,1280) c(n,q)];
% % %     end
% % %     q = find(sort(ds(n,in)-ds(n,1280)) > 0,1,'first');
% % %     if ~isempty(q)
% % %         vp(n,:) = [ds(n,1280) c(n,q)];
% % %     end
end
clear w q in

figure;
imagesc(X); colormap gray; hold on;
j = linspace(0.1*size(X,1),0.9*size(X,1),j);
for n = 1:size(corr,1)
    if vp(n,2)~=0
        plot(cumsum(vp(n,:)),[j(n) j(n)],'-r');
        plot(sum(vp(n,:)),j(n),'>r');
    end
    if vn(n,2)~=0
        plot(cumsum(vn(n,:)),[j(n) j(n)],'-g');
        plot(sum(vn(n,:)),j(n),'<g');
    end
end
