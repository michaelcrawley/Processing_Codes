function schlieren_average()

directory = uigetdir;
flist = dir(directory);
total = zeros(401,1026);
n = 0;

for m = 1:length(flist)
    fprintf([num2str(m) ' of ' num2str(length(flist)) '\n'])
    [path,name,ext] = fileparts([directory filesep flist(m).name]);
    if ~flist(m).isdir&&~isempty(strmatch(ext,'.png'))
        q = double(imread([directory filesep flist(m).name]));
        q = uint8(q/max(q(:))*255);
        total = total + double(q);
        n = n+1;
    end
end
avg = uint8(total/n/max(total(:)/n)*255);

imagesc(avg)
colormap gray
imwrite(avg,[directory filesep 'AVERAGE.png'],'png')