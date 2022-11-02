function [vde, jac] = FokkerPlanckVDE(t, x, h, driftParams, diffusionParams, betamax, calcJacobian)
   % vde = FokkerPlanckVDE(t, u, h, driftParams, diffusionParams, betamax, calcJacobian)
   %
   % RHS of Variational Differential Equations for FokkerPlanckODE
   % containing nominal ODE, state VDE, and parameter VDE
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
   % Nominal ODE:   y' = f(t,y,p)               where   y == u
   % State VDE:    Gy' = df/dy * Gy             where  Gy := dy/dy0
   % Param VDE:    Gp' = df/dy * Gp  +  df/dp   where  Gp := dy/dp
   %
   %
   % The "nominal VDE" with augmented input vector:  uu' = ff(t,uu,p)
   % with:
   %          / u  \                 / f(t,u,p)                   \
   %    uu =  | Gu |    ff(t,uu,p) = | f_u(t,u,p)*Gu              |
   %          \ Gp /                 \ f_u(t,u,p)*Gp + f_p(t,u,p) /
   %
   % Using the method of lines, neglecting time variable t, we have
   %          f =  - a(p) * Au +  D(p) * Bu     - vector 
   %       dfdu =  - a(p) * A  +  D(p) * B      - matrix
   %       dfdp = - da(p) * Au + dD(p) * Bu     - matrix
   %
   % The VDE's jacobian, i.e. the derivative w.r.t. the augmented state uu
   % is 
   %            /  fu(t,u,p)                        0          0      \
   %   dffduu = |  fuu(t,u,p)*Gu                fu(t,u,p)      0      |
   %            \  fuu(t,u,p)*Gp + fpu(t,u,p)       0      fu(t,u,p)  /
   %
   %
   % INPUT:           t --> time
   %                  x --> state vector
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
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
   
   % If calcJacobian is not specified, use default value
   if (nargin < 7), calcJacobian = false; end
   
   % split x into u, Gu, GP parts
   dim = length(x) / 3;
   u  = x(      1:  dim);
   Gu = x(  dim+1:2*dim);
   Gp = x(2*dim+1:3*dim);
   
   % get the required quantities from FokkerPlanckODE:
   if calcJacobian
      [du, dfdu, dfdp, dfdup] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 5);
   else
      [du, dfdu, dfdp]        = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, 4);
   end
   
   % VDE state part:  Gu' = f_u * Gu
   dGu = dfdu * Gu;
   
   % VDE parameter part:  Gp = f_u * Gp + f_p
   dGp = dfdu * Gp + dfdp;
 
   % assemble rhs for full VDE
   vde = [ du ; dGu ; dGp ];
   

   % was the Jacobian requested?
   if (calcJacobian)
      % NOTE: dfduu = 0;
      % TODO: could be beneficial to be set up as sparse in a more clever way
      jac = [ dfdu  ,     0  ,    0 ; ...
                 0  ,  dfdu  ,    0 ; ...
              dfdup ,     0  , dfdu  ];
   end

   
   
   
end