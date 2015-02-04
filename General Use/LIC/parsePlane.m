function [VM,iI,iJ] = parsePlane(V,MS,OL,SW)
%This function takes a plane or set of planes and parses it into pieces
%according to the following:
% 1) "MS" is the maximum size allowed for an individual parse dimension.
% 2) "OL" is the percentage overlap of the parses.
% 3) "SW" is a switch allowing for choice of rectangular or square parses.
%Thus, the program will parse the plane into the minimum number of pieces
%which are smaller than the maximum allowed size with a desired amount of
%overlap. In addition, the user can force the pieces to be square. The
%result is a set of matrices which are all of the same size and will have
%dimensions less than "MS", but also minimizes the number of pieces within
%the constraint of "MS." Since all parses will be the same size, it is
%likely that the total size of the combined matrix will be slightly larger 
%than the input matrix. The overflow will be located on the ends of the
%matrix (e.g. end-2:end).
%
%INPUTS
% V - plane or set of planes. Parsing will occur over the first two
%  dimensions of V and multiple planes can be parsed simultaneously by
%  stacking planes in the third dimension.
% MS - the maximum allowed parse length. NOTE: typical parses will be
%  smaller than "MS."
% OL - The percentage overlap between parses as fraction of parse length.
% SW - Allows user to choose square or rectangular parses. Valid inputs are
%  'square' or 'rect'.
%
%OUTPUTS
% VM - a cell array containing the parsed matrices.
% iI - the indices along the first dimension for the parses. Allows for
%  easy reconstruction of the original matrix. The starting index for each
%  parse is located in the first row of "iI". The ending indices are in the
%  second row.
% iJ - the indices along the second dimension for the parses.


n = size(V,1); n(2) = size(V,2);
if strcmpi(SW,'square') && n(1)~=n(2)
    n = max(n);
    V(n,n,:) = 0;
    n(2) = n;
end

ns = ceil(n/((1-OL)*MS));      %Number of segments in one dimension
nps = ceil(n./(ns*(1-OL)+OL));   %Number of points per segment
nOL = ceil(OL*nps); %Number of points overlap

np = nps + (ns-1).*(nps-nOL);    %Number of points 
V(np(1),np(2),:) = 0;

iI = [1 (nps(1)-nOL(1)+1:nps(1)-nOL(1):np(1)-nOL(1)-1)];    %Starting indices of parses in first dimension
iI(2,:) = iI+nps(1)-1;  %Ending indices of parses in first dimension

iJ = [1 (nps(2)-nOL(2)+1:nps(2)-nOL(2):np(2)-nOL(2)-1)];    %Starting indices of parses in second dimension
iJ(2,:) = iJ+nps(2)-1;  %Ending indices of parses in second dimension

VM = cell(ns);
for j = 1:ns
    for k = 1:ns
        VM{j,k} = V(iI(1,j):iI(2,j),iJ(1,k):iJ(2,k),:);
    end
end
