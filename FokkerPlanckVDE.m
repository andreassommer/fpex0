function dx = FokkerPlanckVDE(t, x, h, driftParams, diffusionParams, betamax)
   % du = FokkerPlanckVDE(t, u, h, driftParams, diffusionParams, betamax)
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
   % INPUT:           t --> time
   %                  x --> state vector
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
   %                        
   %
   % OUTPUT:  df --> vector/matrix containing the requested partial right hand sides
   %                 (nomial ODE, state VDE, parameter VDE)
   %
   %
   % Note: This function is vectorized.
   %
   % If SolvIND is not available, use this function as pure Matlab implementation.
   %
   % Andreas Sommer, Oct2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   % split x into u, Gu, GP parts
   dim = length(x) / 3;
   u  = x(      1:  dim);
   Gu = x(  dim+1:2*dim);
   Gp = x(2*dim+1:3*dim);
   
   % get the required quantities from FokkerPlanckODE:
   [du, dfdu, dfdp] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, true, true, true);

   % VDE state part:  Gu' = f_u * Gu
   dGu = dfdu * Gu;
   
   % VDE parameter part:  Gp = f_u * Gp + f_p
   dGp = dfdu * Gp + dfdp;
 
   % assemble rhs for full VDE
   dx = [ du ; dGu ; dGp ];
   
end