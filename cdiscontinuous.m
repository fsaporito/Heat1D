function y  = cdiscontinuous (x)

% c: oscillation equation -> elastic coefficient
%    heat equation -> thermal conductivity

	% c from 0 to 0.5
	if x < 0.5
	
		y = 1;
	
	% c from 0.5 to 1
	else
	
		y = 10;
		
	end

end
