function out = FokkerPlanckVDE(t, x, h, driftParams, diffusionParams, betamax, np, calcJacobian)
   % vde = FokkerPlanckVDE(t, u, h, driftParams, diffusionParams, betamax, np, calcJacobian)
   %
   % RHS of Parameter-VDE (Variational Differential Equations) 
   % for FokkerPlanckODE containing nominal ODE and parameter VDE
   %
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   %             where a is the driftFcn, and D is the diffusionFcn
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % Jacobian:   df/dpDrift = - d/dpDrift a(t,p) * u_x
   %             df/dpDiff  =   d/dpDiff  D(t,p) * u_xx
   %
   %
   % Nominal ODE:   y' = f(t,y,p)                y(0) = y0(p)     NOTE: y == u
   % State VDE:    Gy' = df/dy * Gy             Gy(0) = I         NOTE: Gy := dy/dy0    - NOT DONE HERE
   % Param VDE:    Gp' = df/dy * Gp  +  df/dp   Gp(0) = dy0/dp    NOTE: Gp := dy/dp
   %
   %
   % The "nominal VDE" with augmented input vector:  uu' = ff(t,uu,p)
   % with:
   %          / u  \                 / f(t,u,p)                   \
   %    uu =  | Gu |    ff(t,uu,p) = | f_u(t,u,p)*Gu              |
   %          \ Gp /                 \ f_u(t,u,p)*Gp + f_p(t,u,p) /
   %
   % Using the method of lines, neglecting time variable t, we have
   %
   %          f =  - a(p) * Au +  D(p) * Bu     - vector 
   %
   %       dfdu =  - a(p) * A  +  D(p) * B      - matrix
   %
   %       dfdp = - da(p) * Au + dD(p) * Bu     - matrix
   %            = - Au * gradient_p[a(p)] + Bu * gradient_p[D(p)]
   %
   %
   % The VDE's jacobian, i.e. the derivative w.r.t. the augmented state uu, formally reads
   % as 
   %            /  fu(t,u,p)                        0          0      \
   %   dffduu = |  fuu(t,u,p)*Gu                fu(t,u,p)      0      |
   %            \  fuu(t,u,p)*Gp + fpu(t,u,p)       0      fu(t,u,p)  /
   % Note: fuu = 0
   % NOTE: THIS IS A TENSOR - NOT A MATRIX
   %
   % 
   %
   %
   % INPUT:           t --> time
   %                  x --> state vector
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
   %                 np --> total number of parameters (FP + IC)
   %       calcJacobian --> flag indicating to compute Jacobian of VDE [default: false]
   %                        
   %
   % OUTPUT:  df --> vector/matrix containing the requested partial right hand sides
   %                 (nomial ODE, state VDE, parameter VDE)
   %
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Andreas Sommer, Oct2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   %
   %
   % NOTES:
   % 
   % Gp    = [ Gp_FP    | Gp_IC    ]
   % Gp(0) = [ dIC/dpFP | dIC/dpIC ] = [ 0 | dIC/dpIC ]    % initial sensitivty only depends on IC-parameters
   %
   % Gp' = f_u * Gp + f_p
   %     = [ df/du ] * [ Gp_FP | Gp_IC ]  +  [ df/dpFP | df/dpIC ] 
   %     = [ df/du ] * [ Gp_FP | Gp_IC ]  +  [ df/dpFP |    0    ]   % no dependency on IC-parameters of rhs f
   %     = [ df/du * Gp_FP + df/dpFP  |  df/du * Gp_IC + 0 ]
   %
   %
   % NOTATION   N = length(u)   n = np = length(p)
   %
   %      / f_1 \             /  df1/dp1  ...  df1/dpn  \
   % f =  |  :  |        fp = |     :             :     |  ==> fp(i,j) = dfi/dpj
   %      \ f_N /             \  dfN/dp1  ...  dfN/dpn  /
   %
   %       N x 1                        N x np
   %
   % Transform fp into a vector:
   % 
   % F(k) = fp(i,j)  with  k = (i-1)*N + j    <-- matlab-like
   % is a vector of length m := N * np
   %
   %          /  dF1/du1  ...  dF1/duN  \
   % dF/du =  |     :             :     |    This is our "jacobian"
   %          \  dFm/du1  ...  dFm/duN  /
   %
   %
   % For the ODE (not the VDE), we have:
   %    f(u,p) = -a(p) * A * u  +  D(p) * B * u
   %   fp(u,p) = - Au * da/dp   +  Bu * dD/dp    ==>  fp(i,j) = - da/dp(i) * (A*u)(j) + dD/dp(i) * (B*u)(j)
   %
   % d/du_k fp(u,p) = - da/dp(i) * d/du (A*u)(j) + dD/dp(i) * d/du (B*u)(j)
   %     fup(i,j,k) = - da/dp(i) * A(j,k)  +  dD/dp(i) * B(j,k)
   %
   
   
   % if calcJacobian is not specified, use default value
   if (nargin < 7), calcJacobian = false; end
   
   % number of parameters
   np_FP = length(driftParams) + length(diffusionParams);
   np_IC = np - np_FP;
   
   % split x into u and Gp parts
   N = length(x) / ( 1 + np );
   u = x(1:N);
   Gp = reshape(x(N+1:end), N, np);
   
   % get the required quantities from FokkerPlanckODE:
   if calcJacobian
      [f, dfdu, dfdpFP, dfdupFP] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 5);
   else
      [f, dfdu, dfdpFP]          = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 4);
   end
   
   % add zeros to build complete dfdp (zeros: no dependency of ODE rhs on IC parameters)
   dfdp = [dfdpFP sparse(N, np_IC)];
   
   % VDE state part:  Gu' = f_u * Gu
   % dGu = dfdu * Gu;
   
   % DEFAULT: Build VDE
   if (~calcJacobian)
      % VDE parameter part:  Gp = f_u * Gp + f_p
      dGp = dfdu * Gp + dfdp;
      dGp = reshape(dGp, np*N, 1);
      % assemble rhs for full VDE
      out = [ f ; dGp ];
      return
   end
   
   % JACOBIAN
   % if we reached this point, calcJacobian must have been true
   
   % was the Jacobian dff/duu == d(VDE)/d(u,Gp) requested?
   if (calcJacobian)
      % add zeros to build complete dfdp (zeros: no dependency of ODE rhs on IC parameters)
      dfdup = [dfdupFP ; sparse(N*N*np_IC, 1)];
      error('Not yet implemented')
      out = 1;
   end

   
   
   
end