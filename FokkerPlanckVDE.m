function out = FokkerPlanckVDE(t, sol_u, Gp, h, driftParams, diffusionParams, betamax, np, calcJacobian)
   % vde = FokkerPlanckVDE(t, sol_u, Gp, h, driftParams, diffusionParams, betamax, np, calcJacobian)
   %
   % RHS of parameter-VDE (Variational Differential Equations) for FokkerPlanckODE.
   % The parameter sensitivity matrix Gp is separated into columns.
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
   % State VDE:    Gy' = df/dy * Gy             Gy(0) = I         NOTE: Gy := dy/dy0
   % Param VDE:    Gp' = df/dy * Gp  +  df/dp   Gp(0) = dy0/dp    NOTE: Gp := dy/dp
   %
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
   % INPUT:           t --> time
   %              sol_u --> nominal solution as matlab sol object
   %                 Gp --> sensitivity matrix as vector (columns)
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
   %                 np --> total number of parameters (FP + IC)
   %       calcJacobian --> flag indicating to compute Jacobian of VDE [default: false] - NOT YET IMPLEMENTED
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
   %         f(u,p) = -a(p) * A * u  +  D(p) * B * u
   %        fp(u,p) = - Au * da/dp   +  Bu * dD/dp    ==>  fp(i,j) = - da/dp(i) * (A*u)(j) + dD/dp(i) * (B*u)(j)
   % d/du_k fp(u,p) = - da/dp(i) * d/du (A*u)(j) + dD/dp(i) * d/du (B*u)(j)
   %     fup(i,j,k) = - da/dp(i) * A(j,k)  +  dD/dp(i) * B(j,k)
   %
   
   
   % if calcJacobian is not specified, use default value
   if (nargin < 7), calcJacobian = false; end
   
   % number of parameters
   np_FP = length(driftParams) + length(diffusionParams);
   np_IC = np - np_FP;
   
   % reshape input vector to Gp matrix for easy matrix-vector multiplications
   N  = length(Gp) / np;  % dimension of the state space
   Gp = reshape(Gp, N, np);
   
   % evaluate nominal solution
   u = deval(sol_u, t);
   
   % get the required quantities from FokkerPlanckODE:
   if calcJacobian
      [dfdu, dfdpFP, dfdupFP] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 7);
   else
      [dfdu, dfdpFP]          = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 6);
   end
   
   % add zeros to build complete dfdp (zeros: no dependency of ODE rhs on IC parameters)
   dfdp = [dfdpFP sparse(N, np_IC)];

   % DEFAULT: Build VDE
   if (~calcJacobian)
      dGp = dfdu * Gp + dfdp;         % VDE parameter part:  Gp = f_u * Gp + f_p
      dGp = reshape(dGp, np*N, 1);    % reshape into vector
      out = dGp;                      % assign rhs output 
      return
   end
   
   % JACOBIAN
   % if we reached this point, calcJacobian must have been true
   
   % was the Jacobian dff/duu == d(VDE)/d(u,Gp) requested?
   if (calcJacobian)
      error('VDE Jacobian computation is not yet implemented.');
      % add zeros to build complete dfdp (zeros: no dependency of ODE rhs on IC parameters)
      dfdup = [dfdupFP ; sparse(N*N*np_IC, 1)];
      out = 1;
      return
   end

   error('Unreachable point reached. What happened?');
   
   
end