function [NOM_sol, SENS_sol] = FPEX0_simulate(FPEX0setup, p)
  % sol = FPEX0_simulate(FPEX0setup, p)
  % [sol, solPSENS] = FPEX0_simumlate(FPEX0setup, p)
  %
  % Simulates Fokker-Planck with specified parameters for 
  % FP drift, FP diffusion, and initial function.
  %
  % INPUT: 
  %         FPEX0setup --> problem setup for FPEX0
  %                  p --> full parameter vector containing parameters for
  %                        Fokker-Planck drift, Fokker-Planck diffusion 
  %                        and initial distribution 
  %
  % OUTPUT: 
  %            sol --> sol object with nominal solution 
  %       solPSENS --> sol object with parameter sensitivities
  %
  % Andreas Sommer, 2017, Aug2022, Nov2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu
  %
  
  % store running index
  persistent runID
 
  % initialize or increment running ID
  if isempty(runID), runID = 0; end
  runID = runID + 1;
  
  % sensitivities requested?
  withSensitivities = false;
  if (nargout >= 2)
     withSensitivities = true;
  end
  
  
  % extract parameters
  p_FPdrift     = FPEX0setup.Parameters.extract_p_FPdrift(p);
  p_FPdiffusion = FPEX0setup.Parameters.extract_p_FPdiffusion(p);
  p_IC          = FPEX0setup.Parameters.extract_p_iniDist(p);
  
  % quick access to temperature grid
  gridT         = FPEX0setup.Grid.gridT;
   
  % retrieve "time" horizon for integrator
  t0tf          = FPEX0setup.Grid.gridTdot([1 end]);  
  
  % generate right hand side, jacobian
  FPrhs         = FPEX0setup.make_FPrhsFcn(p_FPdrift, p_FPdiffusion);
  FPjac         = FPEX0setup.make_FPjacFcn(p_FPdrift, p_FPdiffusion);
    
  % setup integrator and update jacobian therein
  integrator    = FPEX0setup.Integration.integrator;
  odeoptions    = FPEX0setup.Integration.updateODEJacobian(FPjac);
  
  
  % initialize some variables
  NOM_time    = NaN;       SENS_time    = NaN;
  NOM_success = false;     SENS_success = false;
  NOM_sol     = [];        SENS_sol     = [];
  
  
  % make matlab issue an error instead of a warning, if integration fails
  saveWarning = warning('error', 'MATLAB:ode15s:IntegrationTolNotMet'); %#ok<CTPCT> % Undocumented 1st arg set to 'error'
  
  % integrate
  try
     % initial distribution u0 nominal only
     u0 = FPEX0setup.IniDistFcn(gridT, p_IC);
     % measure time for integration
     NOM_time = tic();
     NOM_sol  = integrator(FPrhs, t0tf, u0, odeoptions);
     NOM_time = toc(NOM_time);
     NOM_success = true;
  catch err
     FPEX0_showError('Nominal integration failed!', err);
  end
  
  
  
  % calculate sensitivity
  if withSensitivities
     % Extract things from setup 
     [~, GpIC]     = FPEX0setup.IniDistFcn(gridT, p_IC);    % initial distribution sensitivity GpIC = du0/dp0
     VDEintegrator = FPEX0setup.Integration.VDEintegrator;
     VDErhs = FPEX0setup.make_VDErhsFcn(p_FPdrift, p_FPdiffusion, NOM_sol);  
     N      = FPEX0setup.Grid.N;
     np_FP  = length(FPEX0setup.Parameters.idxFPall);
     % initial condition for VDE
     GpFP = zeros(N, np_FP);     % for FP-parameters, no dependency of u0 on pFP
     Gp0  = [ GpFP , GpIC ];     % initial value GP(0)
     Gp0  = reshape(Gp0, [], 1); % reshape into vector
     % start integration
     try
        SENS_time = tic();
        SENS_sol  = VDEintegrator(VDErhs, t0tf, Gp0, odeoptions);
        SENS_time = toc(SENS_time);
        SENS_success = true;
     catch err
        FPEX0_showError('Sensitivity integration failed!', err);
     end
  end
  
 
  % restore warning condition
  warning(saveWarning); 

  
  
  % store simulation data (for debugging)
  if (FPEX0setup.debugMode.storeSimulation)
     storeID   = sprintf('sim_%d', runID);
     storeData = struct('timeNOMINAL', NOM_time, 'successNOMINAL', NOM_success, 'solNOMINAL', NOM_sol, ...
                        'timeSENS'   , SENS_time   , 'successSENS'   , SENS_success   , 'solSENS'   , SENS_sol    );
     FPEX0setup.store(storeID, storeData);
  end
  
end