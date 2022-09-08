classdef FPEX0_class_setup < handle
% FPEX0_class_setup
%
% Creates FPEX0 configuration.
%
% CONSTRUCTOR:
%
%    setup = FPEX0_class_setup(GridObj, ParametersObj, IntegrationObj, FPdriftFcn, FPdiffusionFcn, IniDistFcn)
%
% INPUT:    GridObj --> grid object of class FPEX0_class_grid
%     ParametersObj --> parameters object of class FPEX0_class_parameters
%    IntegrationObj --> integration settings object of class FPEX0_class_integration
%        FPdriftFcn --> function handle to Fokker-Planck drift
%    FPdiffusionFcn --> function handle to Fokker-Planck diffusion
%        IniDistFcn --> function handle to initial distribution
%
% OUTPUT:   setup --> structure containing configuration
%
% See FPEX0_exampleSetup.m for example usage.
%
% Andreas Sommer, Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


properties (SetAccess = private)
   Grid
   Parameters
   Measurements
   Integration

   Storage = struct();
   
   FPdriftFcn
   FPdiffusionFcn
   IniDistFcn
   
   FPODE_jacPattern
         
end

properties (Access = public)
   debugMode = struct( 'calcresvec' , true );
end


methods
   
   
   % CONSTRUCTOR    
   function setup = FPEX0_class_setup(GridObj, ParametersObj, IntegrationObj, FPdriftFcn, FPdiffusionFcn, IniDistFcn)
      setup.Grid           = GridObj;
      setup.Parameters     = ParametersObj;
      setup.Integration    = IntegrationObj;
      setup.FPdriftFcn     = FPdriftFcn;
      setup.FPdiffusionFcn = FPdiffusionFcn;
      setup.IniDistFcn     = IniDistFcn;
      
      % generate and store the jacobian pattern (does not change as grid does not change)
      setup.FPODE_jacPattern = setup.make_FPODE_jacpattern(setup.Grid.N);
      
      % update the integration options
      opts = odeset( 'JPattern'    , setup.FPODE_jacPattern  ...
                   , 'NonNegative' , 1:setup.Grid.N          );
      setup.Integration.updateOptions(opts);
   end
   
   
   
   % GENERATOR: Right hand side function
   function rhsFcn = make_rhsFcn(obj, p_FPdrift, p_FPdiffusion)
      rhsFcn = @(t,u) FokkerPlanckODE(t, u, obj.Grid.h, obj.FPdriftFcn, p_FPdrift, obj.FPdiffusionFcn, p_FPdiffusion);
   end
   
   
   
   % GENERATOR: Jacobian function
   function jacFcn = make_jacFcn(obj, p_FPdrift, p_FPdiffusion)
      jacFcn = @(t,u) FokkerPlanckODE_dfdu(t, u, obj.Grid.h, obj.FPdriftFcn, p_FPdrift, obj.FPdiffusionFcn, p_FPdiffusion);
   end
   
   
   % Store something
   function store(obj, name, thing)
      % Stores the thing in Storage.(name)
      obj.Storage.(name) = thing;
   end
      
   
   
   % IMPORT MEASUREMENTS
   function importMeasurements(obj, ID, heatrate, T, values, varargin)
      % importMeasurements(FPEX0setup, ID, heatrate, T, values, skip)
      %
      % Imports measurements into FPEX0setup.
      % An interpolation to the integration grid is performed and the
      % interpolated values are stored.
      % The measurement data is used by FPEX0_calcresvec().
      %
      % INPUT: 
      %     FPEX0setup --> FPEX0_class_setup object
      %             ID --> name/identifier of experiment
      %       heatrate --> heating rate
      %              T --> temperature grid / measurement positions
      %         values --> measurement values (e.g. cp-values); 
      %                    where values(k) denotes the measurement value at T(k)
      %       gridskip --> take measurement on every n-th grid point only (default: 1)
      %
      % OUTPUT:    sol --> sol object returned by matlab integrator
      %
      % Andreas Sommer, Aug-Sep2022
      % andreas.sommer@iwr.uni-heidelberg.de
      % code@andreas-sommer.eu
      %

      % check optional skip argument
      gridskip = 1;
      if (nargin >= 5), gridskip = varargin{1}; end

      % ensure that the rate (=time) lies inside the gridTdot
      if (heatrate < obj.Grid.gridTdot(1)) || (heatrate > obj.Grid.gridTdot(end))
         warning('importMeasurements: Specified rate %g outside grid bounds [%g, %g] !', ...
            heatrate, obj.Grid.gridTdot(1), obj.Grid.gridTdot(end));
      end

      % ensure we have a single column
      T      = reshape(T     , [], 1);
      values = reshape(values, [], 1);

      % create a measurements structure
      measstruct.ID        = ID;
      measstruct.heatrate  = heatrate;
      measstruct.rawT      = T;
      measstruct.rawValues = values;

      % interpolate the raw values at the integration grid
      % no extrapolation: set values outside known area to zero
      measstruct.temperatures = obj.Grid.gridT(1:gridskip:end);
      measstruct.values       = interp1(T, values, measstruct.temperatures, 'pchip', 0.0);
      measstruct.values       = reshape(measstruct.values, [] , 1 );

      % add measstruct to measurements
      obj.addMeasStruct(measstruct);

   end
   
   % Add a measstruct
   function addMeasStruct(obj, measstruct)
      obj.Measurements = [obj.Measurements measstruct];
   end
   
   
end




methods (Static)
   
   % HELPER: Jacobian dudf sparsity pattern
   function jacpattern = make_FPODE_jacpattern(N)
      e = ones(N,1);
      jacpattern = spdiags([e e e], -1:1, N, N);
   end

end



% END OF CLASS
end




% setup debug modes
% % FPEX0.debugMode.calcresvec = true;
% % FPEX0.debugMode.simulate   = false;

   
% simulation data
% FPEX0.sim.sol = NaN;   % will be filled by FPEX0_simulate





