classdef sets
	%SETS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
	end
	
	methods (Static)
		[elements, counts] = unique(a);
		[elements, aCount, bCount] = intersect(a, b);
		[c, ia, ib] = union(a, b);
		[elements] = setdiff(a, b);
	end
	
end

