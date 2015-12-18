function [rRaw] = ReadFile(filename,pp)
%Reads in data from file specified in 'filename', based on the information
%supplied in the processing parameters structure, and reshapes as
%necessary.  Outputs raw data (pressure trace vs time for all channels).

%Last updated by Michael Crawley on 2012-05-10

    if pp.isbinary
		fid = fopen([pp.src_dir filesep filename],'r');
		rRaw = fread(fid,'float32'); fclose(fid);
		if pp.TS	%If timestamp is present in files
			rRaw = reshape(rRaw,pp.BS,pp.Nch+1,[]);	%parses data into channels along second dimension and blocks along third dimension
			rRaw = permute(rRaw,[1 3 2]);	%reorders data into channels along third and blocks along second dimension
			rRaw = reshape(rRaw,[],pp.Nch+1);	%reshapes data into 2-D matrix with all blocks for a channel in one column
			rRaw = rRaw(:,2:end);	%removes timestamp
		else
			rRaw = reshape(rRaw,pp.BS,pp.Nch,[]);
			rRaw = permute(rRaw,[1 3 2]);
			rRaw = reshape(rRaw,[],pp.Nch);
		end
	else
		rRaw = dlmread([pp.src_dir filesep filename],pp.delimI); %reads data file
		if pp.TS	%If timestamp is present in files
			rRaw = rRaw(:,2:end);	%removes timestamp
		end
		S = length(rRaw(:,1));
		if round(S/pp.BS)~=S/pp.BS
			error('Data is not an integer number of blocks')
		end
    end
    rRaw = -rRaw; %invert sign of voltage trace. The combination of 4939/2690 mic/conditioner results in a output voltage that is phase-inverted from the original pressure trace.
end