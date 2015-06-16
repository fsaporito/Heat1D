function fhcellarray = trapezoid (N, x, h, fname)

	% Force function definition
	f = str2func(fname);

	% fh Array Definition
	fh = zeros (N-1,1);
	
	% For Loop Trapezoid
	for i=1:N-1
		
		fh(i) = (h(i) + h(i+1))/2 * f(x(i));
	
	end
	
	
	% Method Name
	integralName = 'Trapezoid Method';
	
	
	% Return Data
	fhcellarray = { fh; integralName };

end
