function odecellarray = eulerImplicit (N, x, boundaryCondition, Mh, Kh, fh, dt, Nk, u0name)

	% U0 definition
	u0 = str2func (u0name);
	
	% Matrix that holds all the uhk arrays
	% where uhk(:,k) is the solution at step k (at time k*dt)
	uhk = zeros(N-1,Nk); 

	% uhk(:0) is the initial condition
	% uhk(:,1) must be calculated out of the loop

	% uh0 = uhk(:,0) definition
	uh0 = zeros(N-1,1);

	% uh0 = uhk(:,0) calculation
	for i=1:N-1
		uh0(i) = u0(x(i));
  end
    
  % Neumann Boundary Condition
  if ( strcmp (boundaryCondition, 'DN') == 1)
    uh0 = [ uh0; u0(x(N)) ];
    uhk = zeros(N,Nk); 
  end

	% uh0 = uhk(:,0) ----> uhk(:,1)
	% A is always a constant in our case, but in theory
	% can depend upon the time
	A = (1/dt) * Mh + Kh;
	b = fh + (1/dt) * (Mh * uh0);
	uhk(:,1) = A\b;
    
    
	% uhk(:,k) ----> uhk(:,k+1)
	% A is always a constant in our case,
	% but it can in theory depend upon the time
	for k=1:Nk
		A = (1/dt) * Mh + Kh;
		b = fh + (1/dt) * (Mh * uhk(:,k));
		uhk(:,k+1) = A\b;
	end
	
	ODEmethod = 'Implicit Euler';
	
	odecellarray = { uhk; uh0; ODEmethod };


end
