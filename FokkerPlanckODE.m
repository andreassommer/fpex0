function [out1, out2, out3, out4, out5] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, calcflag)
   % (1) f    = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 1)
   % (2) dfdu = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 2)
   % (3) dfdp = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 3)
   % (4) [du, dfdu, dfdp] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 4)
   % (5) [du, dfdu, dfdp, dfdup] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 5)
   %
   % ODE RHS of Fokker-Planck PDE by using the method of lines (MOL) approach.
   % also including output required for Jacobian computation and VDE computations.
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)  =:  f(t,u,p)
   % 
   % FD-Approx:  Au := u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             Bu := u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %             with matrices A and B being appropriate FD stencils 
   %
   % This function is able to compute several intertwined rhs. Formulation in terms of ODE.
   %
   %     (1)  f = f(t,u,p)
   %          ---  The nominal ODE rhs
   %
   %     (2)  fu = df/du(t,u,p)
   %          ---  The derivative of the ODE rhs w.r.t. state u
   %
   %     (3)  fp = df/dp(t,u,p)
   %          ---  The derivative of the ODE rhs w.r.t. parameter p
   %
   %     (4)  [f, fu, fp] = [f, df/du, df/dp]
   %          ---  Calculates (1)-(3) simultaneously with individual outputs f, dfdu, dfdp
   %
   %     (5)  [f, fu, fp, fup] = [f, df/du, df/dp, df/dup] 
   %          ---  Additionally to (4), also calculates second derivatives df/dup 
   %               NOTE: Second derivative w.r.t. to state u is zero:  f_uu = df/duu = 0
   %
   %
   % INPUT:           t --> time
   %                  u --> state vector
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
   %             calcF --> flag to calculate du    (nominal rhs)                           [default: true ]
   %           calcDFDU --> flag to calculate dfdu  (derivative of rhs w.r.t. states u)     [default: false]
   %           calcDFDP --> flag to calculate dfdp  (derivative of rhs w.r.t. parameters p) [default: false]
   %         calcDVDEDU --> flag to calculate dVDEdu (derivative of VDE w.r.t. state u)     [default: false]
   %
   % OUTPUT:  Call types (1)-(5):
   %                  f --> rhs vector
   %               dfdu --> derivative of rhs w.r.t. state variable u
   %               dfdp --> derivative of rhs w.r.t. drift and diffusion parameters
   %              dfdup --> 2nd order derivative of rhs w.r.t. states and parameters f_up
   %
   %          Call type (5)
   %                VDE --> VDE rhs vector in augmented state uu
   %            dVDEduu --> derivative of VDE rhs w.r.t. the augmented state uu
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Note: - vectorized implementation
   %       - if dfdp or dVDEdu is requested, the internaly used driftFcn and diffusionFcn must be able 
   %         to return their derivatives w.r.t. the parameters as second return argument
   %            df/dpDrift = - d/dpDrift a(t,p) * u_x    % drift
   %            df/dpDiff  =   d/dpDiff  D(t,p) * u_xx   % diffusion
   %
   % Andreas Sommer, Aug2022, Oct2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   % if no calcflag was given, use default value
   if (nargin < 7), calcflag = 1; end

   % initialize flags
   flag__F       = 1;   % f (nominal ODE)
   flag__DFDU    = 2;   % dfdu
   flag__DFDP    = 3;   % dfdp
   flag__DF1_ALL = 4;   % all 1st order derivatives
   flag__DF2_ALL = 5;   % all 1st and nonzero 2nd order derivatives
   
%    flag__VDE_U
%    flag__VDE_P
%    flag__DVDE_DU

   % check what was requested
   switch calcflag
      case flag__F       ,  calcF = true;   calcDFDU = false;  calcDFDP = false;  calcDFDUP = false;  % f (nominal ODE)
      case flag__DFDU    ,  calcF = false;  calcDFDU = true;   calcDFDP = false;  calcDFDUP = false;  % dfdu
      case flag__DFDP    ,  calcF = false;  calcDFDU = false;  calcDFDP = true;   calcDFDUP = false;  % dfdp
      case flag__DF1_ALL ,  calcF = true;   calcDFDU = true;   calcDFDP = true;   calcDFDUP = false;  % f, dfdu, dfdp
      case flag__DF2_ALL ,  calcF = true;   calcDFDU = true;   calcDFDP = true;   calcDFDUP = true;   % f, dfdu, dfdp, dfdup
      otherwise
         msg = 'FokkerPlanckODE: Bad output request.';
         disp(msg);
         error(msg);
   end

   % allocate persistent storage for differentiation stencils
   persistent NN A B
   
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
   
   % vectorized evaluation --- with stored stencil matrices
   for i = 1:vectors
      Au(:,i) = full(A * u(:,i));    % vectors are not sparse anymore
      Bu(:,i) = full(B * u(:,i));    % so change them into full vectors
   end
   
   
   % evaluate drift and diffusion --- possibly with derivates
   if (calcDFDP)
      [alpha, dalpha] = FPEX0_driftFcn(t, driftParams);
      [D    , dD    ] = FPEX0_diffusionFcn(t, diffusionParams, betamax);
      % ensure dalpha and dD are "single-row"-vectors
      dalpha = reshape(dalpha, 1, []);
      dD     = reshape(dD    , 1, []);
   else
      alpha = FPEX0_driftFcn(t, driftParams);
      D     = FPEX0_diffusionFcn(t, diffusionParams, betamax);
   end
   
   
   % calculate what was requested
   if (calcF)    , du    = - alpha * Au +  D * Bu; end   % Vector
   if (calcDFDU) , dfdu  = - alpha * A  +  D * B;  end   % Matrix
   if (calcDFDP) , dfdp  = -dalpha * Au + dD * Bu; end   % Matrix - requires correct shape of dalpha and dD
   if (calcDFDUP), dfdup = -dalpha * A  + dD * B;  end   % Matrix - requires correct shape of dalpha and dD
   
   % assign outputs
   switch calcflag
      case flag__F       ,  out1 = dVDEdu;
      case flag__DFDU    ,  out1 = dfdu;
      case flag__DFDP    ,  out1 = dfdp;
      case flag__DF1_ALL ,  out1 = du;  out2 = dfdu;  out3 = dfdp;
      case flag__DF2_ALL ,  out1 = du;  out2 = dfdu;  out3 = dfdp;  out4 = dfdup;
   end
   
   
end