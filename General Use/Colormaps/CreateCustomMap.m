function [map] = CreateCustomMap(caxis,trange,brange)
    X = 2^10;
    disc = linspace(caxis(1),caxis(2),X);  %discretize color axis range
    [~,I] = min(abs(disc));
    
    %create points for each map section, convert from RGB to 0 to 1
    tmap = makeLine(trange{1},trange{2},I)/255; 
    rmap = makeLine(brange{1},brange{2},X-I)/255;

    %concatenate into one map
    map = [tmap;rmap];    
end