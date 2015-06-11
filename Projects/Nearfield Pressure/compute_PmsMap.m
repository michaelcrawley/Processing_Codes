function [fullmap, submap, supmap, x, y] = compute_PmsMap(src_dir,flist)
%The purpose of this function is to read in filtered, near-field pressure
%data for a given condition, compute the mean-squared pressure, and
%organize the output based on microphone location.

    %Gather data
    for n = 1:length(flist)
        %Read in file
        tmp = load([src_dir filesep flist{n}]);
        if ~isfield(tmp,'matversion') || tmp.matversion < 1.1
            error('Wrong matfile-version');
        end
        
        %Initialize variables during first pass
        if n == 1
            [~,I1,I2] = intersect(roundn(tmp.filter_params.x,-4),roundn(tmp.daq_params.phys.x*tmp.daq_params.phys.D,-4));
            fullmap = zeros(length(flist),length(I1));
            submap = fullmap;
            supmap = fullmap;
            x = fullmap;
            y = fullmap;
        end
        
        %Process relevant data
        x(n,:) = tmp.daq_params.phys.x(I2);
        y(n,:) = tmp.daq_params.phys.y(I2);
        
        intwaveform = reshape(permute(tmp.data.intwaveform,[1 3 2]),[],length(tmp.filter_params.x));
        subsonic = reshape(permute(tmp.data.subsonic,[1 3 2]),[],length(tmp.filter_params.x));
        supersonic = reshape(permute(tmp.data.supersonic,[1 3 2]),[],length(tmp.filter_params.x));
        
        fullmap(n,:) = std(intwaveform(:,I1)).^2;
        submap(n,:) = std(subsonic(:,I1)).^2;
        supmap(n,:) = std(supersonic(:,I1)).^2;        
    end
    
    %Reorganize data
    [~,I1] = sort(y(:,1)); %sort based on array radial position
    y = y(I1,:);
    x = x(I1,:);
    fullmap = fullmap(I1,:);
    submap = submap(I1,:);
    supmap = supmap(I1,:);
end