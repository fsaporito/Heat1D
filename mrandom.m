function meshcellarr = mrandom(N)

% mrandom: generate a random unidimensional mesh
%
% mrandom(N)
%
% N: Nodes' Number

	% Random Mesh (Unique nodes, sorted)
	x = unique(sort(rand(1,N-1)));
	
	% New Nodes'Number (Removed eventual duplicated nodes)
	n = length(x)+1;

	meshtype = 'Random Mesh';

	meshcellarr = { x; n; meshtype};

end
