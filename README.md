# Heat1D
Heat equation solution with finite element method on uniform and random unidimensional mesh

 fem1d: Solve the monodimensional heat equation 
		 rho*u_t - (cu_x)_x = f

		 with Dirichlet Dirichlet Conditions:
		 u(0) = alpha
		 u(1) = beta
		 or
		 with Dirichlet Neumann Conditions:
		 u(0) = alpha
		 c(1)u'(1) = beta
 		 
 Input: 

 fem1d(N, meshname, bctype, bc, fname, cname, rhoname, dt, Tmax, integname, u0name, odename)

 N -> Nodes Number

 meshname -> Function name (without .m) containing the mesh
			  - Uniform mesh, "muniform.m"
			  - Quadratic mesh, "muquadratic.m"
			  - Random mesh, "random.m"

 bctype -> String, 'DD' or 'DN' selects the conditions type

 bc -> Array holding the boundary conditions, two elements
	   - Boundary Condition in 0
	   - Boundary Condition in 1

 fname -> Function name (without .m) containing the f definition 

 cname -> Function name (without .m) containing the c definition 

 rhoname -> Function name (without .m) containing the rho definition 

 dt -> Time step

 Tmax -> Max time

 integname -> Function name (without .m) containing the
			   numerical integration algorithm
			   - Trapezoid method, "trapezoid.m"
			   - Medium point method, "mediumpoint.m"
			   - Simpson Method, "simpson.m"

 u0name -> Function name (without .m) containing the initial data 

 odename -> Function name (without .m) containing the
			 numerical ode solving algorithm
			 - Esplicit Euler, "eulerEsplicit.m"
			 - Implicit Euler, "eulerImplicit.m"


