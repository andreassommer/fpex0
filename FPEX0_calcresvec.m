function [resvec, jacobian] = FPEX0_calcresvec(FPEX0setup, p_all)
   %  [resvec, jacobian] = FPEX0_calcresvec(FPEX0setup, p_all)
   %
   %  INPUT:  
   %     FPEX0setup --> FPEX0 setup class object
   %          p_all --> parameter vector
   %
   %  OUTPUT:   
   %         resvec --> vector of residuals
   %       jacobian --> jacobian matrix (optional, if requested)      [NOT YET IMPLEMENTED]
   %
   %
   %  NOTES: * the settings structure containts the following fields:
   %             .icfun_parameter_count  --> number of parameters of icfun
   %             .gridT                  --> grid for temperatures
   %             .gridTdot               --> grid for temperature gradients
   %         * to avoid implicit weighting, make sure that the computation grid 
   %           and the measurement grid coincide.
   %         * only measurement points that lie on computation grid points are 
   %           used in the calculations
   %
   % Author:  Andreas Sommer, 2017, 2018, 2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   
   
   % jacobian requested?
   if (nargout > 1)
      calcJac = true; 
   else 
      calcJac = false;
   end
   % at the moment...
   if (calcJac == true)
      error('Not yet implemented.');
   end
   
   % debug messages
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
   %        * Ideally, both grids coincide - that avoids implicit weighting.
   
   % get integration "time" grid (heating rates)
   grid_T    = FPEX0setup.Grid.gridT;
   
   % Simlate and store the FP solution, possibly with parameter sensitivites
   [solNominal, solJacobian] = FPEX0_simulate(FPEX0setup, p_all, calcJac);
   
   % evaluate simulation data at heating rates
   simNominal = deval(solNominal, meas_rates);
   
   % start time measurement
   resvecTICid = tic;
     
   % build residual vector
   resvecs = cell(meas_count, 1);
   for k = 1:meas_count
      % select those measurements that lie on a temperature grid point
      [~, compIdxSim, compIdxMeas] = intersect(grid_T, meas_T{k});
      if length(compIdxMeas) ~= length(meas_T{k})
         error('Temperature grid for simulation and measurements do not coincide.');  % for the time being... 
      end
      measVals = meas_values{k}(compIdxMeas);   % measurements restricted to simulation grid
      simVals  = simNominal(compIdxSim, k);        % simulations restricted to measurement grid
      resvecs{k} = measVals - simVals;          % residuals
      if (FPEX0setup.debugMode.showProgress)
         showProgress(meas_T{k}, simVals, measVals, resvecs{k}, meas_count, k);
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




function showProgress(T, simVals, measVals, resVals, n, k)
   persistent hPlots hAxes hFigure
   % create figure if needed;
   if ( isempty(hFigure) || ~isvalid(hFigure) )
      fignum = 339;
      hFigure = figure(fignum);
      hFigure.GraphicsSmoothing = 'off';
   end
   % number n has changed (or first plot): reset
   if (length(hAxes) ~= n)
      hPlots = cell(n,1); 
      hAxes = cell(n,1);
   end
   % Possibly re-use existing plots for performance
   hAx = hAxes{k};
   if ( isempty(hAx) || ~isvalid(hAx) )
      % no axis available, so generate plot
      hAxes{k} = subplot(1, n, k); 
      cla(hAxes{k});          % clear existing content
      hold(hAxes{k}, 'on');   % hold the following plots
      hSimVals  = plot(T, simVals , 'g.-', 'DisplayName', 'Simulation' );
      hMeasVals = plot(T, measVals, 'r.' , 'DisplayName', 'Measurement');
      hResVals  = plot(T, resVals , 'b.' , 'DisplayName', 'Residual'   );
      hPlots{k} = struct('hSimVals', hSimVals', 'hMeasVals', hMeasVals', 'hResVals', hResVals); % store handles
      if (k == n)  % add a legend on last figure
         legend('Location', 'Best');
      end
   else
      % plots available, so just update the data
      hPl = hPlots{k};
      hPl.hSimVals.XData = T;        hPl.hMeasVals.XData = T;         hPl.hResVals.XData = T;
      hPl.hSimVals.YData = simVals;  hPl.hMeasVals.YData = measVals;  hPl.hResVals.YData = resVals;
   end
   % graphics update via drawnow must be done outside!
end


      

