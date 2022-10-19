function sol = FPEX0_simulate(FPEX0setup, p, withDerivatives)
  % sol = FPEX0_simumlate(FPEX0setup, p, withDerivatives)
  %
  % Simulates Fokker-Planck with specified parameters for 
  % FP drift, FP diffusion, and initial function.
  %
  % INPUT: 
  %       FPEX0setup --> problem setup for FPEX0
  %                p --> full parameter vector containing parameters for
  %                      Fokker-Planck drift, Fokker-Planck diffusion 
  %                      and initial distribution 
  %  withDerivatives --> flag indicating to calculate derivatives 
  %                      [default: FALSE]
  %
  % OUTPUT: 
  %         sol --> sol object returned by matlab integrator 
  %
  % Andreas Sommer, 2017, Aug2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu
  %
  
  % store running index
  persistent runID
 
  % initialize or increment running ID
  if isempty(runID), runID = 0; end
  runID = runID + 1;
  
  % extract parameters
  p_FPdrift     = FPEX0setup.Parameters.extract_p_FPdrift(p);
  p_FPdiffusion = FPEX0setup.Parameters.extract_p_FPdiffusion(p);
  p_IC          = FPEX0setup.Parameters.extract_p_iniDist(p);
  
  % evaluate initial distribution
  gridT         = FPEX0setup.Grid.gridT;
  u0            = FPEX0setup.IniDistFcn(gridT, p_IC);
 
  % retrieve "time" horizon for integrator
  t0tf          = FPEX0setup.Grid.gridTdot([1 end]);  
  
  % generate right hand side, jacobian
  FPrhs         = FPEX0setup.make_FPrhsFcn(p_FPdrift, p_FPdiffusion);
  FPjac         = FPEX0setup.make_FPjacFcn(p_FPdrift, p_FPdiffusion);
    
  % setup integrator and update jacobian therein
  integrator    = FPEX0setup.Integration.integrator;
  odeoptions    = FPEX0setup.Integration.updateJacobian(FPjac);
  
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
  
  % store simulation data (for debugging)
  if (FPEX0setup.debugMode.storeSimulation)
     storeID = sprintf('sim_%d', runID);
     storeData = struct('timeSIM', timeSIM, 'success', success, 'sol', sol);
     FPEX0setup.store(storeID, storeData);
  end
  
end