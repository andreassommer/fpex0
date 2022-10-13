function dx = FokkerPlanckVDE(t, x, h, driftFcn, driftParams, diffusionFcn, diffusionParams, nominalODE, stateVDE, paramVDE)
   % dx = FokkerPlanckVDE(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams, what)
   %
   % RHS of Variational Differential Equations for FokkerPlanckODE
   % containing nominal ODE, state VDE, and parameter VDE
   %
   %
   % FP-PDE:     u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
   %             where a is the driftFcn, and b is the diffusionFcn
   % 
   % FD-Approx:  u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
   %             u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
   %
   % Jacobian:  df/dpDrift = - d/dpDrift a(t,p) * u_x      % drift
   %            df/dpDiff  =   d/dpDiff  D(t,p) * u_xx     % diffusion
   %
   %
   % Nominal ODE:                                y' = f(t,y,p)
   % State VDE:    Gy := dy/dy0   fulfills ODE  Gy' = df/dy * Gy
   % Param VDE:    Gp := dy/dp    fulfills ODE  Gp' = df/dy * Gp  +  df/dp
   %
   %
   % INPUT:           t --> time
   %                  x --> state vector
   %                  h --> MOL interval size
   %           driftFcn --> drift function, evaluated at t,driftP
   %        driftParams --> parameter vector for drift function
   %       diffusionFcn --> diffusion function, evaluated at t,driftP
   %    diffusionParams --> parameter vector for diffusion function
   %         nominalODE --> include nominal ODE
   %           stateVDE --> include state VDE
   %           paramVDE --> include parameter VDE
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
   % Andreas Sommer, Sep2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   % persistent storage for efficiency
   persistent NN A B

   
   % separate the input
% % %    compcount = ( nominalODE==true  +  stateVDE==true  +  paramVDE==true);   % ensures binary input interpretation
% % %    dim = length(x) / compcount;
% % %    if compcount >= 1,  u = x(1:dim, :);         end
% % %    if compcount >= 2, Gy = x(dim+1:dim+dim, :); end
% % %    if compcount >= 3, Gp = x(dim+1:dim+dim, :); end
   % NNEEEE !!!!   

   error('Not yet implemented.')
   
   
   % dimension
   N = length(u);
   
   % preparations for concatenation (see bottom)
   dx_ODE      = [];
   dx_stateVDE = []; 
   dx_paramVDE = [];
      

   % evaluate drift and diffusion, derivatives only needed for paramVDE
   if ( paramVDE )
      [a, dadp] = driftFcn(t,driftParams);
      [D, dDdp] = diffusionFcn(t,diffusionParams);
   else
      a = driftFcn(t,driftParams);
      D = diffusionFcn(t,diffusionParams);
   end
   
   
   
   % part shared for nominal ODE and parameter VDE
   % ----------------===========-----=============
   
   if ( nominalODE || paramVDE )
   
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
     
   end
   
   
   
   % part for state VDE and parameter VDE
   % ---------=========-----=============
   if ( stateVDE || paramVDE )
         
      % if dimension has changed, we have to reassemble the stencil matrices
      if (N==NN)
         % Nothing to do. Use existing precomputed matrices.
         % NOTE: Checking N~=NN if NN is empty results always empty array!
         % So we check if they are equal (which is interpreted false if NN is empty),
         % with an empty "true" part and code only in "else" part.
      else
         e1 = ones(N,1);
         e0 = zeros(N,1);
         A = spdiags([e1 ,   e0  , -e1], -1:1, N, N); A(1,2)=0; A(N,N-1)=0;  % 1st order stencil
         B = spdiags([e1 , -2*e1 ,  e1], -1:1, N, N); B(1,2)=2; B(N,N-1)=2;  % 2nd order stencil + robin boundary
         NN = N;
      end
      
      % assemble dfdu
      dfdu = -a * A / (2*h)  +  D * B / h^2;
      
   end
   
   
   
   % part solely for parameter VDE
   % ----------------=============
   if ( paramVDE )

      
   end
   
   
   
   % assemble the individual rhs parts 
   if ( nominalODE )
      dx_ODE = -a * Au + D * Bu;
   end
   if ( stateVDE )
      dx_stateVDE = dfdu;
   end
   if ( paramVDE )
      dx_paramVDE = NaN;
   end
   

   % assemble rhs for full VDE
   dx = [ dx_ODE ; dx_stateVDE ; dx_paramVDE ];
   
   
end