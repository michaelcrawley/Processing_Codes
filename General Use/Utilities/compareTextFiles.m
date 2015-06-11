function dL = compareTextFiles(FN1,FN2)
%Compares two text files line by line and records lines where differences
%occur.

fid = fopen(FN1,'r');
fid2 = fopen(FN2,'r');

n = 1;
d = fgetl(fid);
d2 = fgetl(fid2);

dL = [];
while (sum(d) ~= -1) | (sum(d2) ~= -1)
    if ~strcmp(d,d2)
        dL(end+1) = n;
    end
    n = n+1;
    d = fgetl(fid);
    d2 = fgetl(fid2);
end
fclose(fid);
fclose(fid2);

