classdef parameter
% Describes an experimental parameter.

	properties
		isset = false;
		name = '';
 		symbol = '';
		units = '';
		value = [];
		mod = '';
	end
	
	properties (Dependent)
		print
		describe
	end
	
	methods
		% Constructor
		function p = parameter( name, symbol, units )
			if nargin>0
				p.name = name;
				p.symbol = symbol;
				p.units = units;
			end
		end
		
		% Object set methods
		function obj = set.isset( obj, val )
			if islogical(val),	obj.isset = val;
			else				error('ISSET property must be logical.');
			end
		end
		
		function obj = set.name( obj, val )
			if ischar(val),		obj.name = val;
			else				error('NAME property must be a string.');
			end
		end
		
		function obj = set.symbol( obj, val )
			if ischar(val),		obj.symbol = val;
			else				error('SYMBOL property must be a string.');
			end
		end
		
		function obj = set.units( obj, val )
			if ischar(val),		obj.units = val;
			else				error('UNITS property must be a string.');
			end
		end
		
		function obj = set.value( obj, val )
			if ischar(val) || isnumeric(val)
				obj.value = val;
			else
				error('VALUE property must be numeric or a string.');
			end
			
			obj.isset = true;	% Implicitly true
		end
		
		function obj = set.mod( obj, val )
			if ischar(val),		obj.mod = val;
			else				error('MOD property must be a string.');
			end
		end
		
		% Object get methods
		function val = get.print( obj )
			val = obj.name;
			if ~isempty( obj.symbol )
				val = [val ', ' obj.symbol];
			end
			if ~isempty( obj.value )
				val = [val ': ' num2str(obj.value)];
				if ~isempty( obj.units )
					val = [val obj.units];
				end
			end
		end
		
		function val = get.describe( obj )
			val = obj.name;
			if ~isempty( obj.symbol )
				val = [val ', ' obj.symbol];
			end
			if ~isempty( obj.units )
					val = [val ' [' obj.units ']'];
			end
		end
	end
	
end