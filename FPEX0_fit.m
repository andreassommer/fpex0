function fitsol = FPEX0_fit(FPEX0setup, varargin)
   % fitsol = FPEX0_fit(FPEX0setup, [key-value-pair]* )
   %
   % INPUT:   key-value-pairs:
   %          optimizer --> name of optimizer to be used
   %             fignum --> figure window number to draw in
   %  ICparametrization --> name of parametrization of initial condition (at heat rate 0)
   %
   % OUTPUT: fitsol structure containing solution information
   %
   % Andreas Sommer, 2016-2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   
   % defaults
   optimizer = 'lsqnonlin';

   % process input arguments
   if hasOption(varargin, 'optimizer') ,  optimizer = getOption(varargin, 'optimizer');  end

   
   % integration grid(s)
   gridTdot = FPEX0setup.Grid.gridTdot;
   
   % retrieve parameter values, bounds, indices
   p_0  = FPEX0setup.Parameters.p0;
   p_lb = FPEX0setup.Parameters.p_lb;
   p_ub = FPEX0setup.Parameters.p_ub;
   p_FPdiffIdx = FPEX0setup.Parameters.idxFPdiffusion;

   % Build linear constraints:
   %
   % We want to have non-negative diffusion in whole time horizon
   %     p_FPdiffusion_1 + t * p_FPdfiffusion_2  >=  0   for all t in [0, T]
   %
   % For linear diffusion, it's thus sufficient to enforce that at the grid points 0 and T:
   %     p_FPdiffusion_1 >= 0                            --> via lower bound
   %   - p_FPdiffusion_1 - T * p_FPdfiffusion_2  <=  0   --> as linear constraint
   %
   % Constraint formula for solver fmincon:  A * p <= b
   %
   b_diff_constr = 0;
   A_diff_constr = zeros(1, length(p_0));              % most parameters don't have linear constraints
   A_diff_constr(p_FPdiffIdx) = [-1  -gridTdot(end)];  % only diffusion parameters have linear constraints

   % set function that computes the residual vector
   resvecfun = @(p) FPEX0_calcresvec(FPEX0setup, p);
     
   % optimization
   fprintf('\nStarting:  %s\n', optimizer)
   switch upper(optimizer)
      
      case 'LSQNONLIN'
         % ============= Using SOLVIND derivatives
         lsqnonlin_opts = optimoptions(@lsqnonlin);
         lsqnonlin_opts.SpecifyObjectiveGradient = true;       % jacobian is built by calcresvec
         lsqnonlin_opts.CheckGradients           = false;      % the FD approximation is too poor to verify anything
         lsqnonlin_opts.UseParallel              = true;       % no benefit with analytic jacobian, but it does not hurt
         lsqnonlin_opts.Diagnostics              = 'on';
         lsqnonlin_opts.FunValCheck              = 'off';
         lsqnonlin_opts.StepTolerance            = 1.0d-6;
         lsqnonlin_opts.FunctionTolerance        = 1.0d-10;
         lsqnonlin_opts.OptimalityTolerance      = 1.0;        % 
         lsqnonlin_opts.Display                  = 'iter-detailed';
         lsqnonlin_opts.TypicalX                 = p_0;
         lsqnonlin_opts.SubproblemAlgorithm      = 'factorization';
         lsqnonlin_opts.MaxIterations            = 50;
         lsqnonlin_opts.OutputFcn                = @optimizer_outfun;
         [x,resnorm,res,exit,out,lambda,jac] = lsqnonlin(resvecfun, p_0, p_lb, p_ub, lsqnonlin_opts);
         fitsol = struct('x',x,'resnorm',resnorm,'residual',res,'exitflag',exit,'output',out,...
                         'lambda',lambda,'jacobian',jac,'opts',lsqnonlin_opts);
         displayResult(x, resnorm);

      case 'LSQNONLINFD'
         % ============= Derivative-based solver not a good idea, as FD will probably fail for the PDE solutions
         lsqnonlin_opts = optimoptions(@lsqnonlin);
         lsqnonlin_opts.SpecifyObjectiveGradient = false;      % jacobian is built by calcresvec
         lsqnonlin_opts.CheckGradients           = false;      % the FD approximation is too poor to verify anything
         lsqnonlin_opts.UseParallel              = true;       % very beneficial for FD - initialize pool via parpool(n) 
         lsqnonlin_opts.Diagnostics              = 'on';
         lsqnonlin_opts.FunValCheck              = 'off';
         lsqnonlin_opts.MaxFunctionEvaluations   = 100000;
         lsqnonlin_opts.StepTolerance            = 1.0d-6;
         lsqnonlin_opts.FunctionTolerance        = 1.0d-8;
         lsqnonlin_opts.OptimalityTolerance      = 10.0;       % quite high, but okay for FD
         lsqnonlin_opts.Display                  = 'iter-detailed';
         lsqnonlin_opts.TypicalX                 = p_0;
         lsqnonlin_opts.SubproblemAlgorithm      = 'factorization';
         lsqnonlin_opts.MaxIterations            = 100;        % should take no more than 100 from a good start guess
         lsqnonlin_opts.FiniteDifferenceType     = 'central';  % slower, but accuracy needed
         lsqnonlin_opts.OutputFcn                = @optimizer_outfun;
         [x,resnorm,res,exit,out,lambda,jac] = lsqnonlin(resvecfun, p_0, p_lb, p_ub, lsqnonlin_opts);
         fitsol = struct('x',x,'resnorm',resnorm,'residual',res,'exitflag',exit,'output',out,...
                         'lambda',lambda,'jacobian',jac,'opts',lsqnonlin_opts);
         displayResult(x, resnorm);

         
      case 'FMINCONFD'
         % ============= Using SOLVIND derivatives
         warning('FMINCON is not the method of choice for this kind of problem. Try LSQNONLINFD.')
         %TEST derivatives
         [ff, jj] = scalarobjective(p_all); df = zeros(size(p_all));
         for k=1:length(p_all)
            FDh = 1.0d-7; p_all(k) = p_all(k) + FDh; ff2 = scalarobjective(p_all); 
            p_all(k) = p_all(k) - FDh; df(k) = (ff2 - ff) / FDh;
         end
         keyboard
         fmincon_opts = optimoptions('fmincon');
         fmincon_opts.Algorithm                = 'interior-point';  %'interior-point'; %'sqp'  %'active-set';
         fmincon_opts.SpecifyObjectiveGradient = false;             % no derivatives yet
         fmincon_opts.Display                  = 'iter-detailed';
         fmincon_opts.StepTolerance            = 1.0d-8;
         fmincon_opts.FunctionTolerance        = 1.0d-10;
         fmincon_opts.OptimalityTolerance      = 1.0d-3;            % differs from lsqnonlin
         fmincon_opts.CheckGradients           = true;
         fmincon_opts.TypicalX                 = p_0;
         fmincon_opts.UseParallel              = true;
         fmincon_opts.FiniteDifferenceType     = 'central';         % slower, but accuracy needed
         fmincon_opts.OutputFcn                = @optimizer_outfun;  
         fmincon_opts.MaxIterations            = 1000;
         [x,fval,exit,out,lambda,grad,hessian] ...
            = fmincon(@scalarobjective, p_0, A_diff_constr, b_diff_constr, [], [], p_lb, p_ub, [], fmincon_opts);
         fitsol = struct('x',x,'resnorm',fval,'residual',[],'exitflag',exit,'output',out,...
                         'lambda',lambda,'gradient',grad,'hessian',hessian,'opts',fmincon_opts);
         displayResult(x, fval);
         
      case 'FMINSEARCH'
         % ============ Derivative-free method
         warning('FMINSEARCH cannot ensure non-negativity constraints. Errors to be expected!');
         fminsearch_opts = optimset();
         fminsearch_opts.display     = 'iter';
         fminsearch_opts.TolX        = 1e-4;
         fminsearch_opts.MaxFunEvals = 10000;
         fminsearch_opts.MaxIter     = 10000;
         fminsearch_opts.FunValCheck = 'on';
         fminsearch_opts.UseParallel = 'true';
         fminsearch_opts.Diagnostics = 'on';
         fminsearch_opts.OutputFcn   = @optimizer_outfun;         
         resvecnormfun = @(p) norm(resvecfun(p));
         [x,fval,exit] = fminsearch(resvecnormfun, p_0, fminsearch_opts);
         fitsol = struct('x',x,'fval',fval,'exitflag',exit,'resnorm',fval);
         displayResult(x, fval);
            
      case {'NONE','SIM','SIMULATE'}
         % ============ Simulation
         F = resvecfun(p_0);
         resnorm = norm(F);
         fitsol = struct('x',p_0,'resnorm',resnorm);
         displayResult(p_0, resnorm);
         
      case {'DERIVTEST'}
         % For testing the derivatives:
         [F, J] = resvecfun(p_0);
         for k = 1:length(p_0)
            FD_h = 1.0d-6;                          % FD step size
            FD_p = p_0; FD_p(k) = FD_p(k) + FD_h;   % modified parameter
            FD_F = resvecfun(FD_p);                 % evaluate
            FD_J = (FD_F - F) / FD_h;               % FD approximation 
            figure(22); clf; plot(J(:,k),'-'); hold on; plot(FD_J,'.'); 
            title(sprintf('Displaying derivative for parameter #%d.\n', k));
            drawnow; pause
         end
         resnorm = norm(F);
         fitsol = struct('x',p_0,'resnorm',resnorm);
         displayResult(p_0, resnorm);
         
      otherwise
         error('Unknown optimizer: %s', optimizer);
   end
   
   
   
   % store the solution
   FPEX0setup.store('fitsol',fitsol);
   
   % Confidence regions
   %%% settings.poptci  = FP_calcConfidenceFromSolution(sol, 0.95);

   % show integrals
   %%% FP_checkIntegrals(settings.FPsol, dscdata)
      
   % export (save) results to base workspace and finito
   assignin('base', 'FPEX0_fitsol', fitsol);
   
   % finito
   return
   

   
   % ==================================================
   
   % helpers for fmincon
   function [ff, jj] = scalarobjective(ppp)
      if (nargout==1)
         ff = 0.5 * norm(resvecfun(ppp))^2;
      else
         [FF, JJ] = resvecfun(ppp);
         ff = 0.5 * norm(resvecfun(ppp))^2;
         jj = JJ.' * reshape(FF,[],1);
      end
   end

   

   % output function
   function stop = optimizer_outfun(x,~,state)
      stop = false;
      if strcmp(state,'iter')
         % generate graphic output once (costs: 1 forward simulation)
         FPEX0setup.debugMode.showProgress = true;
         FPEX0_calcresvec(FPEX0setup, x);
         FPEX0setup.debugMode.showProgress = false;
         drawnow('limitrate') % expose/update the graphic
      end
   end


   % display parameter vector
   function displayResult(p,resnorm)
      fprintf('Parameter set:  p = [ ');
      fprintf(' %10.5f ', p);
      fprintf('];   %% residual norm: %g\n', resnorm);
   end


end