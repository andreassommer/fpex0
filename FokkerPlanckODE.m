function [out1, out2, out3, out4] = FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, calcflag)
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
   %            NOTE: df/dup is a tensor. We deliver it as sparse vector dfdup(l) = dfdup(i,j,k)
   %            with l = (i-1)*N^2 + (j-1)*N + k    ( i = 1..np_FP,  j = 1..N,  k = 1..N )
   %            To get the second order mixed derivative d/du d/dp f, observe:
   %                    dfdp(i,j) = -dalpha/dp(i) * Au(j)  + dD/dp(i) * Bu(j)
   %            d/du(k) dfdp(i,j) = -dalpha/dp(i) * A(j,k) + dD/dp(i) * B(j,k)  =: dfdup(i,j,k)
   %            Note that A and B are very sparse, thus also dfdup.
   %          ---  NOTE: Second derivative w.r.t. to state u is zero:  f_uu = df/duu = 0
   %
   %     (6)  [fu, fp] = [df/du, df/dp]
   %          ---  Calculates (2) and (3) simultaneously with individual outputs fu, fp
   %
   %     (7)  [fu, fp, fup] = [df/du, df/dp, df/dup] 
   %          ---  Same as (5) but without (1)
   %
   %
   % INPUT:           t --> time
   %                  u --> state vector
   %                  h --> MOL interval size
   %        driftParams --> parameter vector for drift function
   %    diffusionParams --> parameter vector for diffusion function
   %            betamax --> maximum heat rate (forwarded to FP diffusion function)
   %           calcflag --> flag that determines what to calculate (value range 1 to 5)
   %
   % OUTPUT:  Call types (1)-(5):
   %                  f --> rhs vector
   %               dfdu --> derivative of rhs w.r.t. state variable u
   %               dfdp --> derivative of rhs w.r.t. drift and diffusion parameters
   %              dfdup --> 2nd order mixed derivative information rhs w.r.t. states and parameters, see above
   %
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
   flag__DFDU    = 2;   % df/du
   flag__DFDP    = 3;   % df/dp
   flag__F_DF1   = 4;   % f, df/du, df/dp
   flag__F_DFALL = 5;   % f, df/du, df/dp, (d/du df/dp)
   flag__DF1     = 6;   % df/du, df/dp
   flag__DF2     = 7;   % (d/du df/dp) - mixed second order derivative
   
   

   % check what was requested (grown thing, bitwise flags would have been a better idea)
   switch calcflag
      case flag__F       ,  calcF = true;   calcDFDU = false;  calcDFDP = false;  calcDFDUP = false;  % f (nominal ODE)
      case flag__DFDU    ,  calcF = false;  calcDFDU = true;   calcDFDP = false;  calcDFDUP = false;  % dfdu
      case flag__DFDP    ,  calcF = false;  calcDFDU = false;  calcDFDP = true;   calcDFDUP = false;  % dfdp
      case flag__DF1     ,  calcF = false;  calcDFDU = true;   calcDFDP = true;   calcDFDUP = false;  % dfdu, dfdp
      case flag__DF2     ,  calcF = false;  calcDFDU = false;  calcDFDP = false;  calcDFDUP = true ;  % dfdup
      case flag__F_DF1   ,  calcF = true;   calcDFDU = true;   calcDFDP = true;   calcDFDUP = false;  % f, dfdu, dfdp
      case flag__F_DFALL ,  calcF = true;   calcDFDU = true;   calcDFDP = true;   calcDFDUP = true;   % f, dfdu, dfdp, dfdup
      otherwise
         msg = 'FokkerPlanckODE: Bad output request.';
         disp(msg);
         error(msg);
   end

   % allocate persistent storage for differentiation stencils
   persistent NN A B idxNonzeroA idxNonzeroB nnzA nnzB
   
   % number of grid points and number of simultaneously requested vectors (vectorization)
   [N, vectors] = size(u);
   np_drift     = length(driftParams);
   np_diffusion = length(diffusionParams);
   
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
      idxNonzeroA = find(A);  nnzA = nnz(A);
      idxNonzeroB = find(B);  nnzB = nnz(B);
      NN = N;
   end
   
   % products A*u and B*u are not always needed
   if (calcF || calcDFDP)
      % preallocate vectors for A*u (approx. of u_x) and B*u (approx. of u_xx)
      Au = zeros(N, vectors);
      Bu = zeros(N, vectors);
      % vectorized evaluation --- with stored stencil matrices
      for i = 1:vectors
         Au(:,i) = full(A * u(:,i));    % vectors are not sparse anymore
         Bu(:,i) = full(B * u(:,i));    % so change them into full vectors
      end
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
   if (calcF)   
      % rhs is an Nx1 vector
      f = - alpha * Au +  D * Bu;
   end
   
   if (calcDFDU)
      % state derivative is an NxN matrix
      dfdu  = - alpha * A  +  D * B;
      % NOTE: Actually, 
   end
   
   if (calcDFDP)
      % Matrix - requires correct shape of dalpha and dD
      % dfdp  = - Au * dalpha  + Bu * dD; 
      %   
      dfdp = zeros(N, np_drift + np_diffusion);
      for i = 1:np_drift
         dfdp(:, i) = - dalpha(i) * Au;
      end
      for i = 1:np_diffusion
         dfdp(:, i+np_drift) = dD(i) * Bu;
      end
   end   
      
   
   if (calcDFDUP)
      np_FP = np_drift+np_diffusion;
%       entries = cell(np_FP);
%       indices = cell(np_FP);
%       iii = 1;
%       for i = 1:np_drift
%          entries{iii} = reshape( - dalpha(i) * A , [], 1 );
%          %%%indices{iii} = 1;
%          iii = iii + 1;
%       end
%       for i = 1:np_diffusion
%          entries{iii} = reshape( dD(i) * B, [], 1 );
%          %%%indices{iii} = 1;
%          iii = iii + 1;
%       end
%       
      % Variant 2 -- SLOW
      dfdup = spalloc(N*N*np_FP, 1, nnzA + nnzB);  % preallocate sparse vector
      offset = 0;
      for i = 1:np_drift
         dfdup(offset+idxNonzeroA) = - dalpha(i) * A(idxNonzeroA);
         offset = offset + N*N;
      end
      for i = 1:np_diffusion
         dfdup(offset+idxNonzeroB) = dD(i) * B(idxNonzeroB);
         offset = offset + N*N;
      end
   end
   
   
   % assign outputs
   switch calcflag
      case flag__F       ,  out1 = f;
      case flag__DFDU    ,  out1 = dfdu;
      case flag__DFDP    ,  out1 = dfdp;
      case flag__F_DF1   ,  out1 = f;      out2 = dfdu;  out3 = dfdp;
      case flag__F_DFALL ,  out1 = f;      out2 = dfdu;  out3 = dfdp;  out4 = dfdup;
      case flag__DF1     ,  out1 = dfdu;   out2 = dfdp;
      case flag__DF2     ,  out1 = dfdup;
   end
   
   
end