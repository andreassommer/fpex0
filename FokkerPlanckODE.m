function dx = FokkerPlanckODE(t, u, p, h, driftFcn, diffusionFcn)
   % dx = FokkerPlanckODE(t, u, p, h, driftFcn, diffusionFcn)
   %
   % ODE RHS of Fokker-Planck PDE by using the method of lines (MOL) approach
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % INPUT:       t --> time
   %              x --> state vector
   %              p --> parameter vector for drift and diffusion
   %              h --> MOL interval size
   %       driftFcn --> drift function, evaluated at t,p
   %   diffusionFcn --> diffusion function, evaluated at t,p
   %
   % OUTPUT:  dx --> rhs vector
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Note: Vectorized implementation
   %
   % Andreas Sommer, Aug2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

   % number of grid points and number of simultaneously requested vectors
   [N, vectors] = size(u);

   % preallocate vectors for A*u (approx. of u_x) and B*u (approx. of u_xx)
   Bu = zeros(N, vectors);
   Au = zeros(N, vectors);
   
   % first node  (remember: Matlab is 1-based)
   Bu(1,:) = ( -2*u(1,:) + 2*u(2,:) ) / h^2;
   % Au(1) is zero
   
   % inner nodes (remember: Matlab is 1-based, so starting from 2)
   i = 2:N-1;
   Au(i,:) = ( u(i-1,:) - u(i+1,:) ) / (2*h);           % 1st derivative stencil and scale
   Bu(i,:) = ( u(i-1,:) - 2*u(i,:) + u(i+1,:) ) / h^2;  % 2nd derivative stencil and scale
   
   % last node   (remember: Matlab is 1-based, so last node is N)
   Bu(N,:) = ( -2*u(N,:) + 2*u(N-1,:) ) / h^2;
   % Au(N) is zero
   
   % evaluate drift and diffusion
   alpha = driftFcn(t,p);
   D = diffusionFcn(t,p);
   
   % assemble rhs
   dx = -alpha * Au + D * Bu;   

end