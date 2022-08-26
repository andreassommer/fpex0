function dfdu = FokkerPlanckODE_dfdu(t, u, p, h, driftFcn, diffusionFcn)
   % dfdu = FokkerPlanckODE_dfdu(t, u, p, h, driftFcn, diffusionFcn)
   %
   % Jacobian of FokkerPlanckODE w.r.t. state u, i.e. df/du. 
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % Jacobian:  df/du = -a(t,p) * A + D(t,p) * B
   %            where A is the first-derivative stencil [-1 1]
   %            and B is the second-derivative stencil [1 -2 1] plus Robin boundary
   % 
   % INPUT:       t --> time
   %              x --> state vector
   %              p --> parameter vector for drift and diffusion
   %              h --> MOL interval size
   %       driftFcn --> drift function, evaluated at t,p
   %   diffusionFcn --> diffusion function, evaluated at t,p
   %
   % OUTPUT:  dfdu --> sparse jacobian of FokkerPlanckODE
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Andreas Sommer, Aug2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

   % dimension
   N = length(u);
   
   % evaluate drift and diffusion
   alpha = driftFcn(t,p);
   D = diffusionFcn(t,p);
   
   % assemble Jacobian
   e1 = ones(N,1);
   e0 = zeros(N,1);
   A = spdiags([e1 ,   e0  , -e1], -1:1, N, N); A(1,2)=0; A(N,N-1)=0;
   B = spdiags([e1 , -2*e1 ,  e1], -1:1, N, N); B(1,2)=2; B(N,N-1)=2;

   dfdu = -alpha * A / (2*h)  +  D * B / h^2;

end