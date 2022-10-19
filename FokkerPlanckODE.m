function [out1, out2, out3] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, calcDU, calcDFDU, calcDFDP)
   % du   = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, TRUE, false, false)
   % dfdu = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, false, TRUE, false)
   % dfdp = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, false, false, TRUE)
   % [du, dfdu, dfdp] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, TRUE, TRUE, TRUE)
   %
   % ODE RHS of Fokker-Planck PDE by using the method of lines (MOL) approach.
   % also including output required for Jacobian computation and VDE computations.
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % INPUT:           t --> time
   %                  u --> state vector
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
   %             calcDU --> flag to calculate du    (nominal rhs)                           [default: true ]
   %           calcDFDU --> flag to calculate dfdu  (derivative of rhs w.r.t. states u)     [default: false]
   %           calcDFDP --> flag to calculate dfdp  (derivative of rhs w.r.t. parameters p) [default: false]
   %
   % OUTPUT:    dx --> rhs vector                                              [if calcDU   is true]
   %          dfdu --> derivative of rhs w.r.t. state variable u               [if calcDFDU is true]
   %          dfdp --> derivative of rhs w.r.t. drift and diffusion parameters [if calcDFDP is true]
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Note: - vectorized implementation
   %       - if dfdp is requested, the internaly used driftFcn and diffusionFcn must be able 
   %         to return their derivatives w.r.t. the parameters as second return argument
   %            df/dpDrift = - d/dpDrift a(t,p) * u_x    % drift
   %            df/dpDiff  =   d/dpDiff  D(t,p) * u_xx   % diffusion
   %
   % Andreas Sommer, Aug2022, Oct2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

   % allocate persistent storage for differentiation stencils
   persistent NN A B

   % check optional input args
   if (nargin < 8),  calcDU   = true;  end
   if (nargin < 9),  calcDFDU = false; end
   if (nargin < 10), calcDFDP = false; end

   
   % number of grid points and number of simultaneously requested vectors (vectorization)
   [N, vectors] = size(u);
   
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
      A = A / (2*h);
      B = B / h^2;
      NN = N;
   end
   
   
   % preallocate vectors for A*u (approx. of u_x) and B*u (approx. of u_xx)
   Bu = zeros(N, vectors);
   Au = zeros(N, vectors);
   
   % vectorized evaluation --- direct: Same speed as matrix multiplication below, but we need A and B anyways
   % %    % first node  (remember: Matlab is 1-based)
   % %    Bu(1,:) = ( -2*u(1,:) + 2*u(2,:) ) / h^2;
   % %    % Au(1) is zero
   % %
   % %    % inner nodes (remember: Matlab is 1-based, so starting from 2)
   % %    i = 2:N-1;
   % %    Au(i,:) = ( u(i-1,:) - u(i+1,:) ) / (2*h);           % 1st derivative stencil and scale
   % %    Bu(i,:) = ( u(i-1,:) - 2*u(i,:) + u(i+1,:) ) / h^2;  % 2nd derivative stencil and scale
   % %
   % %    % last node   (remember: Matlab is 1-based, so last node is N)
   % %    Bu(N,:) = ( -2*u(N,:) + 2*u(N-1,:) ) / h^2;
   % %    % Au(N) is zero

   % vectorized evaluation --- with stored stencil matrices
   for i = 1:vectors
      Au(:,i) = A * u(:,i);
      Bu(:,i) = B * u(:,i);
   end
   
   % evaluate drift and diffusion
   if (calcDFDP)
      [alpha, dalpha] = FPEX0_driftFcn(t, driftParams);
      [D    , dD    ] = FPEX0_diffusionFcn(t, diffusionParams, betamax);
   else
      alpha = FPEX0_driftFcn(t, driftParams);
      D     = FPEX0_diffusionFcn(t, diffusionParams, betamax);
   end
   
   % calculate what was requested
   if (calcDU)  , du   = - alpha * Au +  D * Bu; end
   if (calcDFDU), dfdu = - alpha * A  +  D * B;  end
   if (calcDFDP), dfdp = -dalpha * Au + dD * Bu; end
   
   % assign outputs
   if (calcDU && calcDFDU && calcDFDP)
      out1 = du; 
      out2 = dfdu;
      out3 = dfdp;
   elseif (calcDU)
      out1 = du;
   elseif (calcDFDU)
      out1 = dfdu;
   elseif (calcDFDP)
      out1 = dfdp;
   else
      msg = 'FokkerPlanckODE: Bad output request.'; 
      disp(msg);
      error(msg);
   end
   
   
end