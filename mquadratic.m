function meshcellarr = mquadratic(N)

% mquadratic: generate a quadratic unidimensional mesh
%
% mquadratic(N)
%
% N: Nodes' Number

	% X Array Definition
	x = zeros(1,N-1);
	
	% For Loop, x(i) calculation
	for i=1:N-1
		
		x(i) = (i/N)^2;
	
	end
	
	meshtype = 'Quadratic Mesh';

	meshcellarr = { x; N; meshtype};

end
