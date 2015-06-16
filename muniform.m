function meshcellarr = muniform(N)

% mrandom: generate a random unidimensional mesh
%
% mrandom(N)
%
% N: Nodes' Number

	% X Array Definition
	x = zeros(1,N-1);

	% For Loop, x(i) calculation	
	for i=1:N-1
	
		x(i) = i/N;
	
	end
	
	meshtype = 'Uniform Mesh';

	meshcellarr = { x; N; meshtype};

end
