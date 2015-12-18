function corr = xcorr2_1d(X,Y,dim)
% corr = xcorr2_1d(X,Y,dim)
%
%Calculates the normalized 1-dimensional cross correlation of two images
%along the dimension specified by the user. This is determined by summing
%all of the elements in the non-correlated dimension. The sum is then
%normalized by the value at zero displacement.
%
% INPUTS
%X - first image - must be 2-d matrix
%Y - second image - must be same size as X
%dim - dimension of correlation - 1==rows, 2==columns
%
% OUTPUTS
%corr - the correlation profile


if dim==1   %Correlate over the rows
    corr = zeros(size(X,1)*2-1,1);
    for n = 1:size(X,2) %Iterate through all columns
        tmp = xcorr(X(:,n),Y(:,n),'coeff');
        tmp(isnan(tmp)) = 0;
        corr = corr +tmp;
    end
elseif dim==2   %Correlate over the columns
    corr = zeros(1,size(X,2)*2-1);
    for n = 1:size(X,1) %Iterate through all rows
        tmp = xcorr(X(n,:),Y(n,:),'coeff');
        tmp(isnan(tmp)) = 0;
        corr = corr +tmp;
    end
else
    error('Bad dimension')
end

corr = corr/n;
% corr = corr/corr((length(corr)+1)/2);
