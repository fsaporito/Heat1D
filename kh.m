function Kh = kh (N, x, cname, h, bc)

% Kh: Calculate the stiffness matrix with medium point
% C definition
c = str2func(cname);

% Dirichlet Dirichlet Boundary Conditions
if ( strcmp (bc, 'DD') == 1)
	
	% Kh Matrix Definition
	Kh = zeros(N-1,N-1);
	
	% Kh(i) Calculation
	for i=1:N-1
	
		if i==1 % First Row
			xms = ( 0    + x(1))/2;
			xmd = ( x(1) + x(2))/2;
			Kh(1,1) = +c(xms)/h(1) + c(xmd)/h(2);
			Kh(1,2) = -c(xmd)/h(2);
			
		elseif i==N-1 % Last Row (N-1)
			xms = (x(N-2) + x(N-1))/2;
			xmd = (x(N-1) +   1   )/2;
			Kh(N-1,N-2) = -c(xms)/h(N-1);
			Kh(N-1,N-1) = +c(xms)/h(N-1) +c(xmd)/h(N);
			
		else % Generic Row (with 3 non-zero elements)
			xms = ( x(i-1) +  x(i)   )/2;
			xmd = ( x(i)   +  x(i+1) )/2;
			Kh(i,i-1) = -c(xms)/h( i );
			Kh(i,i)   = +c(xms)/h( i ) + c(xmd)/h(i+1);
			Kh(i,i+1) = -c(xmd)/h(i+1);
			
		end
		
	end



% Dirichlet Neumann Conditions
elseif ( strcmp (bc, 'DN') == 1)

	% Kh Matrix Definition
	Kh = zeros(N,N);
	
	% Kh(i) Calculation
	for i=1:N 
	
		if i==1 % First Row
			xms = ( 0    + x(1) )/2;
			xmd = ( x(1) + x(2) )/2;
			Kh(1,1) = +c(xms)/h(1) + c(xmd)/h(2);
			Kh(1,2) = -c(xmd)/h(2);
			
		elseif i==N % Last Row (N)
			xms = ( x(N-1) + x(N) )/2;
			%xmd = ( x(N)   + x(N) )/2;
			Kh(N,N-1) = -c(xms)/h(N);
			Kh(N,N) = +c(xms)/h(N);
						
		else % Generic Row (with 3 non-zero elements)
			xms = ( x(i-1) +  x(i)   )/2;
			xmd = ( x(i)   +  x(i+1) )/2;
			Kh(i,i-1) = -c(xms)/h(i);
			Kh(i,i)   = +c(xms)/h(i) + c(xmd)/h(i+1);
			Kh(i,i+1) = -c(xmd)/h(i+1);
			
		end	
	
	end
	
end
	
	
end
