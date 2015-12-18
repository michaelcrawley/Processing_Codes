function spreadFind()

% q = imread('C:\Documents and Settings\Martin Kearney\My Documents\Schlieren_Heated\20071025\New Folder\AVERAGE.jpg');
[fname,fpath] = uigetfile('*.jpg','Select Image');
q = imread(fullfile(fpath,fname));

S = size(q);
x = [1 S(2)];

column = min(find(q(round(S(1)/2),:)>100))+5;

plot(diff(q(:,column)))
title('Determine Indices for Edge of Flow')
top = input('Enter smaller index: ');
bottom = input('Enter larger index: ');
close

done = false;
imagesc(q)
colormap gray
slope = tan(input('Enter Angle (degrees): ')*pi/180);

while ~done
    py = slope*(x-1) + bottom;
    ny = -slope*(x-1) + top;

    imagesc(q)
    colormap gray
    hold on
    plot(x,py)
    plot(x,ny)
    hold off

    test = input('Enter Angle (degrees) (999 to Escape): ');
    if test == 999
        done = true;
    else
        slope = tan(test*pi/180);
    end
end
fprintf(['The spreading angle is: ' num2str(atan(slope)*180/pi) ' degrees\n'])

