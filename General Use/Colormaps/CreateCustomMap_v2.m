function [map] = CreateCustomMap_v2(cvals,colors,X)
    
    [M,N] = size(colors);
    
    if length(cvals) ~= M, error('Incorrect number of color axis values and correspond color identifiers'); end
    if ~exist('X','var'), X = 255; end

    disc = linspace(cvals(1),cvals(end),X);  %discretize color axis range
    L = cvals(end)-cvals(1);
    
    for n = 1:length(cvals)-1
        %create points for each map section, convert from RGB to 0 to 1
        I = round(X*(cvals(n+1)-cvals(n))/L);
        if I > 0 
            tmap{n} = makeLine(colors(n,:),colors(n+1,:),I)/255;
        else
            tmap{n} = [];
        end
    end

    %concatenate into one map
    map = vertcat(tmap{:});    
end