function [md] = getMetaFromStr(str,setup,md)
% Returns the variables of the format [tag][num][mod] from a string. The 
% string must not contain a path or file extension.

if fileparts(str), error('Input cannot contain path.'); end
if nargin<3, md = struct(); end

% Get reference metadata structure
ref = initParams(setup);
refnames = fieldnames(ref);

%%-----------------------------------------------------------------------%%

% Extract substrings divided by underscores
ss = regexpi( str, '_', 'split' );

% Interpret substrings
for i=1:length(ss)
	% See 'Regular Expressions', 'Tokens', and 'Named Capture' to 
	% understand whats going on here.
	
	% Pull variables of the format '[tag][num][mod]'
	expr = ['^' ...							% Beginning of line
		'(?<tag>[a-zA-Z]+)?' ...			% Match upper and lower
		'(?<num>[0-9\().-]+)?' ...				% Match numerals and decimals
		'(?<mod>[a-zA-Z]+)?' ...			% Match upper and lower
		'$'];								% End of line
	a = regexpi( ss{i}, expr, 'names' );
	
	if ~isempty( a.tag )
		if isfield( md, a.tag )
			var = a.tag;
		else
			if isfield( ref, a.tag )
				% Use an exact reference match if possible
				var = a.tag;
				md.(var) = ref.(var);
			else
				% Otherwise look for an inexact match
				match = strcmpi( a.tag, refnames );
				
				if any(match)
					% Use an inexact reference match if one is found
					var = refnames{match};
					md.(var) = ref.(var);
				else
					% Otherwise create an empty parameter
					var = a.tag;
					md.(var) = parameter;
					warning( 'Unrecognized parameter: %s', var );
				end
			end	
		end
		
		% Explicity declare the variable as being set
		md.(var).isset = true;
	else
		% Skip to next substring if no valid tag is found
		warning( 'Invalid or missing variable tag: %s', ss{i} );
		continue;
	end
	
	if ~isempty( a.num )
		% Set the value
        if isnan(str2double(a.num))
            md.(var).value = a.num;
        else
            md.(var).value = str2double( a.num );
        end

		% Set the modifier if one exists
		if ~isempty( a.mod )
			md.(var).mod = a.mod;
		end
	end
end