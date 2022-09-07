function [resvec, jacobian] = FPEX0_calcresvec(FPEX0setup, p_all)
   %  [resvec, jacobian] = FPEX0_calcresvec(FPEX0setup, p_all)
   %
   %  INPUT:  
   %     FPEX0setup --> FPEX0 setup class object
   %          p_all --> parameter vector
   %
   %  OUTPUT:   
   %         resvec --> vector of residuals
   %       jacobian --> jacobian matrix (optional, if requested)  % NOT AVAILABLE IN MATLAB-ONLY VERSION
   %
   %
   %  NOTE:  the settings structure containts the following fields:
   %            .icfun_parameter_count  --> number of parameters of icfun
   %            .gridT                  --> grid for temperatures
   %            .gridTdot               --> grid for temperature gradients
   %
   %
   % Author:  Andreas Sommer, 2017, 2018, Aug2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   
   
   % sorry, no derivatives available in Matlab-only implementation yet
   if (nargout > 1), error('Derivatives not yet available.'); end
   
   % Debug messages
   if FPEX0setup.debugMode.calcresvec
      paramString = sprintf('%22.16f ', p_all );
      fprintf('calcresvec:  p = [ %s ]\n', paramString);
   end

   
   % For higher heating rates we have less data (sample time points) than for low rates.
   % To avoid implicit weighting, use the gridT and interpolate/evaluate the dsc-data at these grid points!

   % get the measurement data to compare the simulation with
   meas_count  = length(FPEX0setup.Measurements);
   meas_values = { FPEX0setup.Measurements.values       }; % cell array of double arrays
   meas_T      = { FPEX0setup.Measurements.temperatures }; % cell array of double arrays
   meas_rates  = [ FPEX0setup.Measurements.heatrate     ]; % double array
   % NOTES: * The measurement temperatures must be "compatible" with the integration grid,
   %          i.e. they must be a subset of the integration space (=temperature) grid points
   %        * Same holds for measurement rates: they must be a subset of the 
   %          integration time grid.
   
   % get integration "time" grid (heating rates)
   grid_T    = FPEX0setup.Grid.gridT;
   
   % simulate and store the FP solution
   sol = FPEX0_simulate(FPEX0setup, p_all);
   
   % extract simulation data at heating rates
   simdata = deval(sol, meas_rates);
   
   % start time measurement
   resvecTICid = tic;
     
   % build residual vectors
   resvecs = cell(meas_count, 1);
   for k = 1:meas_count
      % select those measurements that lie on a temperature grid point
      [~, compIdxSim, compIdxMeas] = intersect(grid_T, meas_T{k});
      if length(compIdxMeas) ~= length(meas_T{k})
         error('Grid does not fit');  % for the time being... 
      end
      measVals = meas_values{k}(compIdxMeas);   % measurements restricted to simulation grid
      simVals  = simdata(compIdxSim, k);        % simulations restricted to measurement grid
      resvecs{k} = measVals - simVals;          % residuals
      if FPEX0setup.debugMode.calcresvec
         showDEBUG(meas_T{k}, simVals, measVals, resvecs{k}, meas_count, k); %% DEBUG: illustration
      end
   end

   
   % concatenate to residual vector
   resvec = vertcat( resvecs{:} );

   % stop time measurement
   time_resvec = toc(resvecTICid);
   
   % store some statistics in settings closure
   % FPEX0.sim.time_resvec = time_resvec;
   
   % finito
   return

end



function showDEBUG(T, simVals, measVals, resVals, n, k)
   fignum = 339;
   figure(fignum);
   subplot(1, n, k); plot(T, simVals, 'g.-', T, measVals, 'r.', T, resVals, 'b.');
end


      