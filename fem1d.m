function [] = fem1d(N, meshname, bctype, bc, fname, cname, rhoname, dt, Tmax, integname, u0name, odename)

% fem1d: Solve the monodimensional heat equation 
%		 rho*u_t - (cu_x)_x = f
%
%		 with Dirichlet Dirichlet Conditions:
%		 u(0) = alpha
%		 u(1) = beta
%		 or
%		 with Dirichlet Neumann Conditions:
%		 u(0) = alpha
%		 c(1)u'(1) = beta
% 		 
% Input: 
%
% fem1d(N, meshname, bctype, bc, fname, cname, rhoname, dt, Tmax, integname, u0name, odename)
%
% N -> Nodes Number
%
% meshname -> Function name (without .m) containing the mesh
%			  - Uniform mesh, "muniform.m"
%			  - Quadratic mesh, "muquadratic.m"
%			  - Random mesh, "random.m"
%
% bctype -> String, 'DD' or 'DN' selects the conditions type
%
% bc -> Array holding the boundary conditions, two elements
%	   - Boundary Condition in 0
%	   - Boundary Condition in 1
%
% fname -> Function name (without .m) containing the f definition 
%
% cname -> Function name (without .m) containing the c definition 
%
% rhoname -> Function name (without .m) containing the rho definition 
%
% dt -> Time step
%
% Tmax -> Max time
%
% integname -> Function name (without .m) containing the
%			   numerical integration algorithm
%			   - Trapezoid method, "trapezoid.m"
%			   - Medium point method, "mediumpoint.m"
%			   - Simpson Method, "simpson.m"
%
% u0name -> Function name (without .m) containing the initial data 
%
% odename -> Function name (without .m) containing the
%			 numerical ode solving algorithm
%			 - Esplicit Euler, "eulerEsplicit.m"
%			 - Implicit Euler, "eulerImplicit.m"
%



% Mesh definition
mesh = str2func(meshname);

% Force function definition
f = str2func(fname);

% C definition
c = str2func(cname);

% Rho definition
rho = str2func(rhoname);

% Numerical integration method
integral = str2func (integname);

% U0(t) definition
% u0 = str2func (u0name);

% Numerical ODE solutiom method
ode = str2func (odename);

% Boundary Conditions
alpha = 0;
beta = 0;
gamma = 0;

alpha = bc(1);

if (strcmp (bctype, 'DD') == 1) % DD Dirichlet Dirichlet
	beta = bc(2);
end

if (strcmp (bctype, 'DN') == 1) % DN Dirichlet Neumann
	gamma = bc(2);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X Definition

meshcellarr = mesh (N); % Mesh Calculation

x = meshcellarr {1, 1}; % Assign The Calculated x Array
N = meshcellarr {2, 1}; % Updating N
meshname = meshcellarr {3, 1}; % Mesh Type

% If Dirichlet Neumann Add Last Node
if ( strcmp (bctype, 'DN') == 1)
	
	x = [x 1];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h Step Array
% h(i) = x(i) - x(i-1)

h = zeros(1,N); % h Array Definition

% For Loop, h Calculation
for i=1:N
    if i==1 % h(1) = x(1)-x(0), x(0) = 0
		h(1) = x(1)-0;
    elseif i==N % h(N) = x(N)-x(N-1), x(N) = 1
		h(N) = 1-x(N-1);        
    else
		h(i) = x(i)-x(i-1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kh Matrix

Kh = kh (N, x, cname, h, bctype);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fh Array

fhcellarr = integral (N, x, h, fname);

fh = fhcellarr {1, 1};
integmethod = fhcellarr {2, 1};

% Non Homogenous Dirichlet Dirichlet Boundary Conditions
if ( strcmp (bctype, 'DD') == 1)

	% Edit Fh For Non Omogenous Dirichlet Conditions In 0
	xms=( 0   + x(1))/2;
	fh(1) = fh(1) - (-alpha/h(1))*c(xms);
	
	% Edit Fh For Non Omogenous Dirichlet Conditions In 1
	xmd=( x(N-1) + 1 )/2;
	fh(N-1) = fh(N-1) - (-beta/h(N))*c(xmd);
	
end

% Non Homogenous Dirichlet Neumann Boundary Conditions
if ( strcmp (bctype, 'DN') == 1)

	% Edit Fh For Non Omogenous Dirichlet Conditions In 0
	xms=( 0   + x(1))/2;
	fh(1) = fh(1) - (-alpha/h(1))*c(xms);
	
	% Edit Fh For Neumann Conditions In 1
	% fhlast = fh(N) = h(N)/2 * f(x(N), Add Element To Fh Array
	fhlast = ( h(N)/2 ) * f(x(N));
	fh = [ fh; fhlast ];
	
	% Edit Fh For Non Omogenous Neumann Conditions In 1
	% fh(N) = fh(N) + gamma
	fh(N) = fh(N) + gamma;
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear System Solution

uh = Kh\fh;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution Plot

figure(1);

% Plot
lw=2; % LineWidth
if ( strcmp (bctype, 'DD') == 1) % DD plot
	plot([0 x 1],[alpha uh' beta],'.-b','LineWidth',lw)
	hold on
	title ([integmethod,' - DD']);
	xlabel(['x axis - ', meshname]);
	ylabel('y axis - Uh(x)');
	legend('Stationary Solution');
end

if ( strcmp (bctype, 'DN') == 1) % DN plot
	plot([0 x],[alpha uh'],'.-b','LineWidth',lw)
	hold on;
	title ([integmethod,' - DN']);
	xlabel(['x axis - ', meshname]);
	ylabel('y axis - Uh(x)');
	legend('Stationary Solution');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mh Matrix

Mh = zeros(N-1,N-1); % Mh Matrix Definition
	
% For Loop, Mh(i) Calculation
% Use Trapezoid method to obtain a diagonal matrix
for i=1:N-1
	Mh(i,i) = (h(i)+h(i+1))/2 * rho(x(i));
end
	
if ( strcmp (bctype, 'DN') == 1)
	
	Mh = zeros(N,N); % Mh Matrix Definition
	
	% For Loop, Mh(i) Calculation
	% Use Trapezoid method to obtain a diagonal matrix
	for i=1:N-1
		Mh(i,i) = ( h(i) + h(i+1) )/2 * rho(x(i));
	end
	
	% Last Element Calculation
	Mh(N,N) = h(N)/2 * rho(x(N));
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Solution

Nk = round(Tmax/dt); % Step number (use round to obtain an integer)

odecellarr = ode (N, x, bctype, Mh, Kh, fh, dt, Nk, u0name);

uhk = odecellarr {1,1}; 
uh0 = odecellarr {2,1}; 
ODEmethod = odecellarr {3,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Solution: Bidimensional Plot

figure(2);

% LineWidth
lw = 2;

if ( strcmp (bctype, 'DD') == 1)

	% Initial data plot
	plot([0 x 1],[alpha uh0' beta],'.-k','LineWidth',lw)
	hold on;

	% Stationary solution plot
	plot([0 x 1],[alpha uh' beta],'.-b','LineWidth',(lw+5))

	% uhk plot
	for k=1:Nk
		plot([0 x 1],[alpha uhk(:,k)' beta],'.-r','LineWidth',lw)
	end	

	% Title with ODE method name
	title (strcat (ODEmethod,' - DD'));

	% Axis
	xlabel(['x axis - ', meshname]);
	ylabel('y axis - Uh(x)');

	% Legend
	legend ('Initial Data','Stationary Solution','Uhk');

end

if ( strcmp (bctype, 'DN') == 1)
    
    % Initial data plot
    plot([0 x],[alpha uh0'],'.-k','LineWidth',lw)
    hold on;
    
    % Stationary solution plot
    plot([0 x],[alpha uh'],'.-b','LineWidth',(lw+5))

    % uhk plot
    for k=1:Nk
        plot([0 x],[alpha uhk(:,k)'],'.-r','LineWidth',lw)
    end

   

    % Title with ODE method name
	  title (strcat (ODEmethod,' - DN'));

	  % Axis
	  xlabel(['x axis - ', meshname]);
	  ylabel('y axis - Uh(x)');

	  % Legend
	  legend ('Initial Data','Stationary Solution','Uhk');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Solution: Tridimensional Plot


% Tridimensional Plot, Dirichlet Dirichet Condition
if ( strcmp (bctype, 'DD') == 1)

	figure(3);
	
	% General Case (i Generic)
	for k=0:Nk
		for i=1:N-2
			
			X = [x(i)      x(i+1)       x(i+1)       x(i)];
			T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
        
			if (k == 0) % Time = 0
				U = [uh0(i)    uh0(i+1)     uhk(i+1,k+1) uhk(i,k+1)];
			else % General Time Case
				U = [uhk(i,k)  uhk(i+1,k)   uhk(i+1,k+1) uhk(i,k+1)];
			end
        
			patch(X,T,U,'w')
    
		end    
	end

	% i = 0 x_0 = 0
	i = 0;
	for k=0:Nk        
    
		X = [0         x(i+1)       x(i+1)       0];
		T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
    
		if (k == 0) % Time = 0
			U = [alpha    uh0(i+1)     uhk(i+1,k+1) alpha];
		else % General Time Case
			U = [alpha    uhk(i+1,k)   uhk(i+1,k+1) alpha];
		end
    
		patch(X,T,U,'w')

	end

	% i = N-1 x_N = 1
	i = N-1;
	for k=0:Nk        
    
		X = [x(i)      1             1            x(i)];
		T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
    
		if (k == 0) % Time = 0
			U = [uh0(i)    beta     beta    uhk(i,k+1)];
		else % General Time Case
			U = [uhk(i,k)  beta   beta    uhk(i,k+1)];
		end
    
		patch(X,T,U,'w')
   
	end
	
	% Force Tridimensional View
	view(3)

	% Axis
	xlabel(['x axis - ', meshname]);
	ylabel('t - time');
	zlabel('y axis - Uhk(x)');

	% Title with ODE method name
	title ([ ODEmethod, ' - DD', '  Patch View' ]);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Surf Plot
  
  figure(4);
  
  % Spatial coordinates vector, adds 0 and 1 points to node vector
  X = [0 x 1];
  
  % Time coordinates vector definition
  T = zeros (Nk+1,1);
  
  % Time coordinates vector calculation (Uniform time mesh)
  for i=1:Nk+1
     
     T(i) = dt*(i-1);
  
  end
  
  % Uh function matrix definition
  U = zeros(Nk+1,N+1);
  
  % Uh function matrix calculation
  for t=1:Nk+1 % Time for loop
  
    for s=1:N+1 % Space for loop
    
      if (s-1 == 0) % Apply Dirichlet boundary condition in s = 0
        U(t,s) = alpha;
        
      elseif ( s == N+1 ) % Apply Dirichlet boundary condition in s = 1
        U(t,s) = beta;
     
      else
        
        if ( t-1 == 0 ) % Use initial data in t = 0
          U(t,s) = uh0(s-1);
        
        else % General calculation
          U(t,s) = uhk(s-1,t-1);
        
        end
      
      end
      
    end
  
  end
  
  
  % Surface drawing
  surf(X,T,U)
  
  % Coloring
  colormap(jet(N))
  colorbar
  
  % Force Tridimensional View
	view(3)

	% Axis
	xlabel(['x axis - ', meshname]);
	ylabel('t - time');
	zlabel('y axis - Uhk(x)');

	% Title with ODE method name
	title ([ ODEmethod, ' - DD', '  Surface View' ]);
  
	
end

% Tridimensional Plot, Dirichlet Neumann Condition
if ( strcmp (bctype, 'DN') == 1)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Patch Plot
  
  figure(3);
  
	% General Case (i Generic)
	for k=0:Nk
		for i=1:N-1
			
			X = [x(i)      x(i+1)       x(i+1)       x(i)];
			T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
        
			if (k == 0) % Time = 0
				U = [uh0(i)    uh0(i+1)     uhk(i+1,k+1) uhk(i,k+1)];
			else % General Time Case
				U = [uhk(i,k)  uhk(i+1,k)   uhk(i+1,k+1) uhk(i,k+1)];
			end
           
	    patch(X,T,U,'w')
    
	  end    
	end

	% i = 0 x_0 = 0
	i = 0;
	for k=0:Nk        
    
		X = [0         x(i+1)       x(i+1)       0];
		T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
    
		if (k == 0) % Time = 0
			U = [alpha    uh0(i+1)     uhk(i+1,k+1) alpha];
		else % General Time Case
			U = [alpha    uhk(i+1,k)   uhk(i+1,k+1) alpha];
		end
    
    patch(X,T,U,'w')

	end

	% i = N x_N = 1
	i = N;
	for k=0:Nk        
    
		X = [x(i)      1             1            x(i)];
		T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
    
		if (k == 0) % Time = 0
			U = [uh0(i)    uh0(i)     uhk(i,k+1)    uhk(i,k+1)];
		else % General Time Case
			U = [uhk(i,k)  uhk(i,k)   uhk(i,k+1)    uhk(i,k+1)];
		end
    
    patch(X,T,U,'w')
   
	end
  
  
  % Force Tridimensional View
	view(3)


	% Axis
	xlabel(['x axis - ', meshname]);
	ylabel('t - time');
	zlabel('y axis - Uhk(x)');

	% Title with ODE method name
	title ([ ODEmethod, ' - DN', '  Patch View' ]);
	figure(3);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Surf Plot
  
  figure(4);
  
  % Spatial coordinates vector, adds 0 point to node vector
  X = [0 x];
  
  % Time coordinates vector definition
  T = zeros (Nk+1,1);
  
  % Time coordinates vector calculation (Uniform time mesh)
  for (i=1:Nk+1)
     
     T(i) = dt*(i-1);
  
  end
  
  % Uh function matrix definition
  U = zeros(Nk+1,N+1);
  
  % Uh function matrix calculation
  for (t=1:Nk+1) % Time for loop
  
    for (s=1:N+1) % Space for loop
    
      if (s-1==0) % Apply Dirichlet boundary condition in s = 0
        U(t,s) = alpha;
     
      else
        
        if (t-1==0) % Use initial data in t = 0
          U(t,s) = uh0(s-1);
        
        else % General calculation
          U(t,s) = uhk(s-1,t-1);
        
        end
      
      end
      
    end
  
  end
  
  
  % Surface drawing
  surf(X,T,U)
  
  % Coloring
  colormap(jet(N))
  colorbar
  
  % Force Tridimensional View
	view(3)

	% Axis
	xlabel(['x axis - ', meshname]);
	ylabel('t - time');
	zlabel('y axis - Uhk(x)');

	% Title with ODE method name
	title ([ ODEmethod, ' - DN', '  Surface View' ]);
	

	
	
	
end


end % EndFunction
