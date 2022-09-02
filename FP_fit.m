function results = FP_fit(varargin)
   % results = FP_fit( [key-value-pair]* )
   %
   % INPUT:   key-value-pairs:
   %          optimizer --> name of optimizer to be used
   %             fignum --> figure window number to draw in
   %  ICparametrization --> name of parametrization of initial condition (at heat rate 0)
   %
   % OUTPUT: results structure
   %

   % access global configuration
   global FPEX0
   
   % defaults
   optimizer = 'lsqnonlinFD';
   iniDist   = 'frasersuzuki';
% % %    draw_sol_fignum   = 666;
% % %    FPwithDrain       = false;   

   % process input arguments
   if hasOption(varargin, 'optimizer') ,   optimizer = getOption(varargin, 'Optimizer');  end
   if hasOption(varargin, 'iniDist')   ,     iniDist = getOption(varargin, 'iniDist');    end
% % %    if hasOption(varargin, 'fignum')    , draw_sol_fignum = getOption(varargin, 'fignum');            end
     
   % integration grid(s)
   gridTdot = FPEX0.grid.gridTdot;
% % %    gridStep = 0.1;
% % %    gridT    = 70:gridStep:180;
% % %    gridTdot = [0 ; reshape( sort(unique(arrayfun(@(x) x.Tinfo.Tstep, dscdata))), [], 1) ];
% % %    %gridT    = reshape(gridT, [], 1);
   
   
   % !!! WICHTIG !!!
   % Die Diffusion und Drift (ihr Zahlenwert) haengt mit der Gitterweite zusammen!
   % Deswegen alle Simulationen immer bei identischem gridT machen!
   
   
   
   
   switch lower(iniDist)
      case {'frasersuzuki',''}
         % parameters of initial distribution, and its bounds  FRASER-SUZUKI
         %         r  ;   h ;  z  ;   wr  ;   sr
         p0    = [ 2  ;  70 ; 132 ;   3   ;  0.5  ];
         p0_lb = [ 1  ;  15 ; 110 ;   0.1 ;  1d-5 ];
         p0_ub = [ 3  ; 150 ; 150 ;  50.0 ;  1.00 ];
         icfun = @frasersuzuki;
      case {'hvl','haarhoffvdlinde','haarhoffvonderlinde'} 
         % parameters of initial distribution, and its bounds  HAARHOFF-VON DER LINDE
         %  a0 (peak area), a1 (center of Gaussian part), a2 (SD of Gaussian part), a3 (peak distortion)
         p0    = [ 228  ;  110  ;  3 ;   -1.5 ];
         p0_lb = [ 160  ;   30  ;  1  ;  -2.0 ];
         p0_ub = [ 250  ;  200  ; 20  ;   2.0 ];
         icfun = @haarhoffvdlinde;
      otherwise
         error('Unknown ICparametrization: "%s"', iniDist);
   end

   
   % parameters for Fokker-Planck
   %        [    p_FPdrift       ;      p_FPdiff    ];
   pFP    = [   0.1   ;  +0.3    ;  0.1   ;  0.5    ];  % lineare drift und diffusion
   pFP_lb = [  -0.1   ;  -0.1    ;  0.00  ;  0.00   ];  %
   pFP_ub = [  10.1   ;  10.1    ;  1.01  ;  1.01   ];  %

   % TEST DEBUG: Konstante Diffusion
   % pFP_lb(4) = 0; pFP_ub(4) = 0; 
   
   % Beachte: Finale Diffusion als Parameter!  
   % Formel: Diffusion = pD1 + t*(pD2-pD1)/betamax   
   %   wobei  betamax: maximale heizrate / zeit
   %          pD1, pD2: erster und zweiter p_FPdiff parameter
   % 
   % -> Lineare Diffusion plus Bedingung, dass Diffusion am Ende nicht-negativ ist.
   
   % TEST DEBUG: Konstante Diffusion und Drift
   % p_lb([2 4]) = 0; p_ub([2 4]) = 0; 
  
   % TEST DEBUG: Mindestdiffusion
   % p_lb(3) = 0.000001;   
   
% % %    % (re-)scale Fokker-Planck-Parameters
% % %    scaleFPp = 10000;
% % %    pFP    = scaleFPp * pFP;
% % %    pFP_lb = scaleFPp * pFP_lb;
% % %    pFP_ub = scaleFPp * pFP_ub;
   
   % all parameters
   p_all    = [pFP    ; p0   ];
   p_all_lb = [pFP_lb ; p0_lb];
   p_all_ub = [pFP_ub ; p0_ub];
  

   % Linear constraints for Diffusion
   % Diffusion parameters are on position 3 & 4 in parameter vector
   %     p_FP_diffusion_1 + t * p_FP_dfiffusion_2  >=  0   for all t in [0, T]
   % So it's sufficient to ensure that at the grid points 0 and T:
   %     p_FP_diffusion_1 >= 0                             --> via lower bound
   %   - p_FP_diffusion_1 - T * p_FP_dfiffusion_2  <=  0   --> as linear constraint
   A_diff_constr      = zeros(1, length(p_all));
   A_diff_constr(3:4) = [-1  -gridTdot(end)];    % diffusion parameters at idx 3 & 4
   b_diff_constr      = 0;

%    % setup settings
%    settings.np0               = np0;
%    settings.np                = np;
%    settings.p0_lb             = p0_lb;
%    settings.p0_ub             = p0_ub;
%    settings.p_lb              = p_lb;
%    settings.p_ub              = p_ub;
%    
%    settings.icfun             = icfun;
%    
%    settings.gridT             = gridT;
%    settings.gridTdot          = gridTdot;
%    settings.FPwithDrain       = FPwithDrain;
%    settings.dscdata           = dscdata;
%    
%    settings.time_FPsim        = [];
%    settings.time_resvec       = []; 
%    settings.time_drawsol      = []; 
%    settings.draw_sol_fignum   = draw_sol_fignum;
%    settings.needNewIntegrator = true;
%    
%    % wrap in closure
%    settingsClosure = makeClosure(settings);
%    

   % set residual function
   resvecfun = @(p) FPEX0_calcresvec(p);
   
     
   % optimization
   fprintf('\nStarting:  %s\n', optimizer)
   switch upper(optimizer)
      
      case 'LSQNONLINFD'
         % ============= Derivative-based solver not a good idea, as FD will probably fail for the PDE solutions
         lsqnonlin_opts = optimoptions(@lsqnonlin);
         lsqnonlin_opts.SpecifyObjectiveGradient = false;
         lsqnonlin_opts.CheckGradients           = false;
         lsqnonlin_opts.UseParallel              = false;
         lsqnonlin_opts.Diagnostics              = 'on';
         lsqnonlin_opts.FunValCheck              = 'on';
         lsqnonlin_opts.MaxFunctionEvaluations   = 100000;
         lsqnonlin_opts.StepTolerance            = 1.0d-6;
         lsqnonlin_opts.FunctionTolerance        = 1.0d-10;
         lsqnonlin_opts.OptimalityTolerance      = 1.0d-1;
         lsqnonlin_opts.Display                  = 'iter-detailed';
         lsqnonlin_opts.TypicalX                 = p_all;
         lsqnonlin_opts.SubproblemAlgorithm      = 'factorization';
         lsqnonlin_opts.MaxIterations            = 1000;
         lsqnonlin_opts.FiniteDifferenceType     = 'central';
         lsqnonlin_opts.OutputFcn                = @optimizer_outfun;
         [x,resnorm,res,exit,out,lambda,jac] = lsqnonlin(resvecfun, p_all, p_all_lb, p_all_ub, lsqnonlin_opts);
         sol = struct('x',x,'resnorm',resnorm,'residual',res,'exitflag',exit,'output',out,'lambda',lambda,'jacobian',jac,'opts',lsqnonlin_opts);
         %%%exportToBaseWorkspace(sol,'FP_fit_sol');
         
      case 'LSQNONLIN'
         % ============= Using SOLVIND derivatives
         lsqnonlin_opts = optimoptions(@lsqnonlin);
         lsqnonlin_opts.SpecifyObjectiveGradient = true;
         lsqnonlin_opts.CheckGradients           = false;
         lsqnonlin_opts.UseParallel              = false;
         lsqnonlin_opts.StepTolerance            = 1.0d-6;
         lsqnonlin_opts.FunctionTolerance        = 1.0d-10;
         lsqnonlin_opts.OptimalityTolerance      = 1.0d-1;
         lsqnonlin_opts.Display                  = 'iter-detailed';
         lsqnonlin_opts.TypicalX                 = p_all;
         lsqnonlin_opts.SubproblemAlgorithm      = 'factorization';
         lsqnonlin_opts.MaxIterations            = 50;
         [x,resnorm,res,exit,out,lambda,jac] = lsqnonlin(resvecfun, p_all, p_all_lb, p_all_ub, lsqnonlin_opts);
         sol = struct('x',x,'resnorm',resnorm,'residual',res,'exitflag',exit,'output',out,'lambda',lambda,'jacobian',jac,'opts',lsqnonlin_opts);
         %%%exportToBaseWorkspace(sol,'FP_fit_sol');
         
      case 'FMINCON'
         % ============= Using SOLVIND derivatives
         warning('This is not working well...')
         % TEST derivatives
         % [ff, jj] = scalarobjective(p_all); df = zeros(size(p_all));
         % for k=1:length(p_all)
         %    FDh = 1.0d-7; p_all(k) = p_all(k) + FDh; ff2 = scalarobjective(p_all); p_all(k) = p_all(k) - FDh; df(k) = (ff2 - ff) / FDh;
         % end
         % keyboard
         fmincon_opts = optimoptions('fmincon');
         fmincon_opts.Algorithm                = 'sqp';%'interior-point'; %'sqp'  %'active-set';
         fmincon_opts.SpecifyObjectiveGradient = true;
         fmincon_opts.Display                  = 'iter-detailed';
         fmincon_opts.StepTolerance            = 1.0d-6;
         fmincon_opts.FunctionTolerance        = 1.0d-10;
         fmincon_opts.OptimalityTolerance      = 1.0d-2; 
         %fmincon_opts.CheckGradients           = true;
         fmincon_opts.TypicalX                 = p_all;
         fmincon_opts.UseParallel              = true;
         [x,fval,exit,out,lambda,grad,hessian] = fmincon(@scalarobjective, p_all, A_diff_constr, b_diff_constr, [], [], p_all_lb, p_all_ub, [], fmincon_opts);
         %[x,fval,exit,out,lambda,grad,hessian] = fmincon(@scalarobjective, p_all, [], [], [], [], [], [], [], fmincon_opts);
         sol = struct('x',x,'resnorm',fval,'residual',[],'exitflag',exit,'output',out,'lambda',lambda,'gradient',grad,'hessian',hessian,'opts',fmincon_opts);
         %%%exportToBaseWorkspace(sol,'FP_fit_sol');
         
         
      case 'FMINSEARCH'
         % ============ Derivative-free method
         fminsearch_opts = optimset('display','iter','TolX',1e-3,'MaxFunEvals',10000,'MaxIter',10000,'FunValCheck','on'); % options
         resvecnormfun = @(p) norm(resvecfun(p));
         [x,fval,exit] = fminsearch(resvecnormfun, p_all, fminsearch_opts);
         sol = struct('x',x,'fval',fval,'exitflag',exit,'resnorm',fval);
         %%%exportToBaseWorkspace(sol,'FP_fit_sol');
            
      case {'NONE','SIM','SIMULATE'}
         % ============ Simulation
         F = resvecfun(p_all);
         resnorm = norm(F);
         sol = struct('x',p_all,'resnorm',resnorm);
         fprintf('Parameter set:  p = [ ');
         fprintf(' %10.5f ', p_all);
         fprintf('];   %% residual norm: %g\n', resnorm);
         
      case {'DERIVTEST'}
         % For testing the derivative:
         [F, J] = resvecfun(p_all);
         for k = 1:length(p_all)
            FD_h = 1.0d-6; FD_p = p_all; FD_p(k) = FD_p(k) + FD_h; FD_F = resvecfun(FD_p); FD_J = (FD_F - F) / FD_h;
            figure(22); clf; plot(J(:,k),'-'); hold on; plot(FD_J,'.'); drawnow;
            fprintf('Displaying derivative for parameter #%d.\n', k); 
            pause
         end
         resnorm = norm(F);
         sol = struct('x',p_all,'resnorm',resnorm);
         fprintf('Parameter set:  p = [ ');
         fprintf(' %10.5f ', p_all);
         fprintf('];   %% residual norm: %g\n', resnorm);
      otherwise
         error('Unknown optimizer: %s', optimizer);
   end
   
   
   
   % Confidence regions
   %%% settings.poptci  = FP_calcConfidenceFromSolution(sol, 0.95);

   % show integrals
   %%% FP_checkIntegrals(settings.FPsol, dscdata)
      
   % export (save) results to base workspace and finito
   %%%results = settings;
   %%%assignin('base', 'FP_fit_result', results);
   
   % finito
   return
   

   
   % ==================================================
   
   % helper for fmincon
   function [ff, jj] = scalarobjective(ppp)
      if (nargout==1)
         ff = 0.5 * norm(resvecfun(ppp))^2;
      else
         [F, J] = resvecfun(ppp);
         ff = 0.5 * norm(resvecfun(ppp))^2;
         jj = J.' * reshape(F,[],1);
      end
   end

   
end



% output function
function stop = optimizer_outfun(x,optimValues,state)
   stop = false;
   if strcmp(state,'iter')
      global FPEX0
      FPEX0.debugMode.calcresvec = true;
      FPEX0_calcresvec(x); % generates graphic output
      FPEX0.debugMode.calcresvec = false;
   end
end

