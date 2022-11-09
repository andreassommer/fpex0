function measstruct = FPEX0_importMeasurements(ID, heatrate, T, values, varargin)
  % measstruct = FPEX0_importMeasurements(ID, heatrate, T, values, skip)
  %
  % Imports measurements into global FPEX0 configuration. 
  % An interpolation to the integration grid is performed and the 
  % interpolated values are stored.
  % The measurement data is used by FPEX0_calcresvec().
  %
  % INPUT:    ID --> name/identifier of experiment
  %     heatrate --> heating rate
  %            T --> temperature grid / measurement positions
  %       values --> measurement values (e.g. cp-values); where
  %                  values(k) denotes the measurement value at T(k)
  %         skip --> take measurement on every n-th grid point only (default: 1)
  % OUTPUT:  sol --> sol object returned by matlab integrator 
  %
  % Andreas Sommer, Aug2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu
  %

% access global configuration
global FPEX0

% check optional skip argument
skip = 1;
if (nargin >= 5), skip = varargin{1}; end

% extract integration space (=temperature) grid and time (=heatrate) grid
gridT    = FPEX0.grid.gridT(1:skip:end);
gridTdot = FPEX0.grid.gridTdot;

% ensure that the rate (=time) lies inside the gridTdot
if (heatrate < gridTdot(1)) || (heatrate > gridTdot(end))
   warning('FPEX0_importMeasurements: Specified rate %g outside grid bounds [%g, %g] !', ...
           heatrate, gridTdot(1), gridTdot(end))
end

% ensure we have a single column
T      = reshape(T     , [], 1);
values = reshape(values, [], 1);

% add measurements
measstruct.ID        = ID;
measstruct.heatrate  = heatrate;
measstruct.rawT      = T;
measstruct.rawValues = values;

% interpolate the raw values at the integration grid
% no extrapolation: set values outside known area to zero
measstruct.temperatures = gridT;
measstruct.values = reshape( interp1(T, values, gridT, 'pchip', 0.0) , [] , 1 );

% add measstruct to FPEX0 global data
FPEX0.measurements = [FPEX0.measurements measstruct];


end