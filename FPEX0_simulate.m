function sol = FPEX0_simulate(p)
  % sol = FPEX0_simumlate(p)
  %
  % Simulates Fokker-Planck with specified parameters for FP drift, diffusion, and initial function.
  %
  % INPUT:     p --> full parameter vector for Fokker-Planck drift and diffusion and initial distribution 
  %
  % OUTPUT:  sol --> sol object returned by matlab integrator 
  %
  % Andreas Sommer, 2017, Aug2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu
  %
  
  % store running index
  persistent runID
  
  % access global configuration
  global FPEX0

  % initialize or increment running ID
  if isempty(runID), runID = 0; end
  runID = runID + 1;
  
  % extract parameters
  p_FPdrift     = FPEX0.parameters.extract.FPdrift(p);
  p_FPdiffusion = FPEX0.parameters.extract.FPdiffusion(p);
  p_IC          = FPEX0.parameters.extract.iniDist(p);
  
  % evaluate initial distribution
  gridT         = FPEX0.grid.gridT;
  u0            = FPEX0.functions.iniDist(gridT, p_IC);
 
  % retrieve "time" horizon for integrator
  t0tf          = FPEX0.grid.gridTdot([1 end]);  
  
  % generate right hand side, jacobian
  driftFcn      = FPEX0.functions.drift;
  diffusionFcn  = FPEX0.functions.diffusion;
  h             = FPEX0.grid.h;
  FPrhs         = FPEX0.functions.make_rhs(h, driftFcn, p_FPdrift, diffusionFcn, p_FPdiffusion);
  FPjac         = FPEX0.functions.make_jac(h, driftFcn, p_FPdrift, diffusionFcn, p_FPdiffusion);
    
  % setup integrator and update jacobian therein
  integrator    = FPEX0.integration.integrator;
  odeoptions    = FPEX0.integration.updateOptionsJacobian(FPEX0.integration.options, FPjac);
  
  % make matlab issue an error instead of a warning, if integration fails
  saveWarning = warning('error', 'MATLAB:ode15s:IntegrationTolNotMet'); %#ok<CTPCT> % Undocumented 1st arg set to 'error'
  
  % start integrattion
  success = false;
  try
     timeSIM = tic();
     sol = integrator(FPrhs, t0tf, u0, odeoptions);
     timeSIM = toc(timeSIM);
     success = true;
  catch err
     disp('Integration failed!');
     fprintf('Error id    : %s\n', err.identifier)
     fprintf('Error string: %s\n', err.message)
     timeSIM = NaN;
     sol = NaN;
  end
    
  % Restore warning condition
  warning(saveWarning);
  
  % store simulation data (kind of DEBUG)
  FPEX0.store(runID, timeSIM, success, sol);
  
end