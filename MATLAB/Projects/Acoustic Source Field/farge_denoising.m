

qmf = MakeONFilter('Battle',5);
N = numel(test);
maxitr = N;

wc = FWT_PO(test,1,qmf);
sigma = sum(abs(wc(:)).^2)/N; %note that this has already been squared
thresh = sqrt(2*sigma*log2(N));
chk = abs(wc) <= thresh;
num = sum(chk(:));

for n = 1:maxitr
    wc_I = zeros(size(wc));
    wc_I(chk) = wc(chk);
    sigma = sum(abs(wc_I(:)).^2)/N;
    thresh = sqrt(2*sigma*log(N));
    chk = abs(wc_I) <= thresh;
    
    if sum(chk(:)) == num
        break;
    else
        num = sum(chk(:));
    end
end

denoised = IWT_PO(wc_I,1,qmf);

% pcolor(denoised);shading interp;colormap jet