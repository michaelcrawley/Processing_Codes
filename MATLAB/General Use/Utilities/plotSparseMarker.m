function plotSparseMarker(x,y,mrk,Nm,labels)
% Creates plots with occasional markers - assumes linear grid
% x - horizontal axis - program can handle a single set of x-data, an array
%   of x-data (one column for each y column) or a cell array (one cell entry
%   per y). Length of x need not be the same for all data sets
% y - array of vertical axis arranged (line1,line2,...) where lineN is a
% column vector of same length as x or cell array.
% mrk - the array of marker choice and color - can handle N-by-2 character
%   array or N-by-2 cell array containing color selection including any color
%   triplet.
%          b     blue          .     point              
%          g     green         o     circle          
%          r     red           x     x-mark             
%          c     cyan          +     plus                 
%          m     magenta       *     star            
%          y     yellow        s     square
%          k     black         d     diamond
%                              v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram
%
% Example: x = (1:10); y = [(1:10)',(11:20)']; mrk = ['bo';'rx'];
%
% Nm - the number of markers on the curve
% labels - the name of each curve for the legend as cell array of strings
%

% colors = 'bgrcmyk';
% markers = '.ox+*sdv^<>ph';
% 
% Narg = nargin;
% if Narg==1
%     y = varargin{1};
%     S = size(y);
%     x = (0:S(1));
%     if S(2)/6 > length(markers)
%         error('Too many data arrays')
%     end
%     for n = 1:S(2)
%         m(n,:) = [color(mod(n,6)+1) markers(ceil(n/6))];
%     end
% elseif Narg==2
%     for n = 1:Narg
%         t(n) = ischar(varargin{n});
%     end
%     if sum(t) > 0
%         x = (0:S(1));
%     else
%         
%     end
% end

    %Reorganizes color/marker combinations to conform with function standards
colors = 'bgrcmyk';
markers = '.ox+*sdv^<>ph';
if ~iscell(mrk)
    q = strfind(colors,mrk(1,1));
    if ~isempty(q)
        C = 1;
        NC = 2;
    else
        C = 2;
        NC = 1;
    end
    S = size(mrk);
    for n = 1:S(1)
        mc{n} = mrk(n,C);
        ms{n} = mrk(n,NC);
    end
    clear mrk; mrk = [ms(:) mc(:)];
else
    if ~ischar(mrk{1,1})
        mrk = fliplr(mrk);
    else
        q = strfind(markers,mrk{1,1});
        if isempty(q)
            mrk = fliplr(mrk);
        end
    end
end

    %Reorganizes data to conform with function standards - assumes each
    %data array is a column in y if y is a matrix
Sy = size(y);
if ~iscell(y)
    Sx = size(x);
    if (Sx(2) > 1)&(Sx(1) > 1)
        for n = 1:Sx(2)
            xt{n} = x(:,n);
            yt{n} = y(:,n);
        end
    else
        for n = 1:Sy(2)
            yt{n} = y(:,n);
            xt{n} = x;
        end
    end
    x = xt; y = yt; clear xt yt;
end

S = size(y);
gcf;
hold on
for n = 1:S(2)
    Rx = max(x{n}) - min(x{n});   %x-range of data
    dx = mean(diff(x{n})); %resolution of data
    L = Rx/Nm;  %x-distance between markers
    N = ceil(L/dx);   %number of points between markers
    Lc = (N:N:length(x{n})-N);
    Lc = round(Lc + (rand(size(Lc))-1/2)*N*0.8);
    Lc = Lc(Lc > 0); Lc = Lc(Lc <= length(x{n}));
    Lc = [1 Lc length(x{n})];
    plot(x{n}(Lc),y{n}(Lc),'LineStyle','none',...
        'Marker',mrk{n,1},'MarkerEdgeColor',mrk{n,2},'MarkerFaceColor',mrk{n,2})
end
legend(labels)
for n = 1:S(2)
    plot(x{n},y{n},'Color',mrk{n,2})
end
hold off


