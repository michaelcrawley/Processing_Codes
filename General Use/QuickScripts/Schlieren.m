L = [100 1125 486-200 486+200];

prefix = 'F:\Research\Mach1.3\Schlieren\20091114\1615\Cropped\M13_T155_i';
suffix = '.png';
N = 50;

for n = 1:50
    d = imread([prefix num2str(n) suffix]);
    imwrite(d(L(3):L(4),L(1):L(2)),[prefix num2str(n) suffix],'png')
end



prefix = 'F:\Research\Mach1.3\Schlieren\20091120\TTR15\F070\m1_ph000\M13_m1_F070_T171_ph000_i';
suffix = '.png';
N = 50;
T = zeros(960,1280);
for n = 1:50
    d = imread([prefix num2str(n) suffix]);
    T = T +double(d);
end
T = uint8(T/N);
clear suffix N

CDF = cumsum(imhist(d)); CDF = CDF/CDF(end);
U = min(find(CDF > 0.999));
dn = uint8(double(d)*255/U);

CDF = cumsum(imhist(T)); CDF = CDF/CDF(end);
U = min(find(CDF > 0.999));
dn = uint8(double(T)*255/U);

dn2 = imrotate(dn,2);
dnT = dn2(200:800,1:1100);
imshow(dnT)

imwrite(dnT,'AVERAGE_m1_F7.5_T171_PH000.png','png')
imwrite(dnT,'INSTANT_m1_F7.5_T171_PH000.png','png')

