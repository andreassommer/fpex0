function dfdp = FokkerPlanckODE_dfdp(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams)
   % dfdu = FokkerPlanckODE_dfdu(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams)
   %
   % Jacobian of FokkerPlanckODE w.r.t. drift and diffusion parameters, i.e. df/dp. 
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   %             where a is the driftFcn, and b is the diffusionFcn
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % Jacobian:  df/dpDrift = - d/dpDrift a(t,p) * u_x    % drift
   %            df/dpDiff  =   d/dpDiff  D(t,p) * u_xx   % diffusion
   % 
   % INPUT:           t --> time
   %                  x --> state vector
   %                  h --> MOL interval size
   %           driftFcn --> drift function, evaluated at t,driftP
   %        driftParams --> parameter vector for drift function
   %       diffusionFcn --> diffusion function, evaluated at t,driftP
   %    diffusionParams --> parameter vector for diffusion function
   %
   % OUTPUT:  dfdp --> sparse jacobian of FokkerPlanckODE
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Andreas Sommer, Sep2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   error('Not yet implemented.')
   
   % dimension
   N = length(u);
   
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
   
   
  
   % evaluate drift and diffusion with derivatives
   [~, dadp] = driftFcn(t,driftParams);
   [~, dDdp] = diffusionFcn(t,diffusionParams);
      
   % assemble rhs for parameter VDE
   dx = -a * Au + D * Bu;   
   
end