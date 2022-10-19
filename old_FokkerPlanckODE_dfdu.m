function dfdu = FokkerPlanckODE_dfdu(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams)
   % dfdu = FokkerPlanckODE_dfdu(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams)
   %
   % Jacobian of FokkerPlanckODE w.r.t. state u, i.e. df/du. 
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   %             where a is the driftFcn, and b is the diffusionFcn
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % Jacobian:  df/du = -a(t,p) * A + D(t,p) * B
   %            where A is the first-derivative stencil [-1 0 1]
   %            and B is the second-derivative stencil [1 -2 1] plus Robin boundary
   % 
   % INPUT:           t --> time
   %                  x --> state vector
   %                  h --> MOL interval size
   %           driftFcn --> drift function, evaluated at t,driftP
   %        driftParams --> parameter vector for drift function
   %       diffusionFcn --> diffusion function, evaluated at t,driftP
   %    diffusionParams --> parameter vector for diffusion function
   %
   % OUTPUT:  dfdu --> sparse jacobian of FokkerPlanckODE
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Andreas Sommer, Aug2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

   persistent NN A B
   
   % dimension
   N = length(u);
   
   % if dimension has changed, we have to reassemble the stencil matrices
   if (N==NN)
      % Nothing to do. NOTE: Checking N~=NN if NN is empty results always empty array!
      % So we check if they are equal (which is interpreted false if NN is empty), 
      % with an empty "true" part and code only in "else" part.
   else
      e1 = ones(N,1);
      e0 = zeros(N,1);
      A = spdiags([e1 ,   e0  , -e1], -1:1, N, N); A(1,2)=0; A(N,N-1)=0;  % 1st order stencil
      B = spdiags([e1 , -2*e1 ,  e1], -1:1, N, N); B(1,2)=2; B(N,N-1)=2;  % 2nd order stencil + robin boundary
      NN = N;
   end
   
   % evaluate drift and diffusion
   %a = driftFcn(t, driftParams);
   %D = diffusionFcn(t, diffusionParams);
   betamax = 20; %%%DEBUG
   a = FPEX0_driftFcn(t, driftParams);
   D = FPEX0_diffusionFcn(t, diffusionParams, betamax);
   
   % assemble the jacobian
   dfdu = -a * A / (2*h)  +  D * B / h^2;

end