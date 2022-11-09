classdef FPEX0_class_integration < handle
% setup = FPEX0_class_integration(varargin)
%
% Generates integrator setup for FPEX0.
%
% CONSTRUCTOR:
%
% setup = FPEX0_class_integration()
% setup = FPEX0_class_integration(integrator)
% setup = FPEX0_class_integration(integrator, options)
% setup = FPEX0_class_integration(integrator, options, VDEintegrator)
% setup = FPEX0_class_integration(integrator, options, VDEintegrator, VDEoptions)
%
% INPUT:  integrator --> handle to Matlab integrator    (default: ode15s)
%            options --> integrator options from odeset (default: see source)
%      VDEintegrator --> handle to Matlab integrator    (default: ode23)
%         VDEoptions --> integrator options from odeset (default: see source)
%
% OUTPUT: setup --> integator setup
%
% Author:  Andreas Sommer, Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%   
   

properties
   
   % ODE: 
   % The initial condition usually has a sharp slope, such that
   % the problem becomes stiff and an implicit solver is needed
   integrator = @ode15s;   
   options    = odeset( 'AbsTol'      , 1e-12    ...
                      , 'RelTol'      , 1e-6     ...
                      , 'BDF'         , 'off'    ...
                      , 'Vectorized'  , 'on'     ...  % ODE rhs is vectorized
                      , 'stats'       , 'off'  );

   % VDE:
   % The sensitivities are usually smooth, and an explicit integrator
   % is the method of choice. 
   % Then, no jacobian information is required.
   VDEintegrator = @ode45;
   VDEoptions = odeset( 'AbsTol'      , 1e-8     ...
                      , 'RelTol'      , 1e-8     ...  % high accuracy for derivatives
                      , 'Vectorized'  , 'off'    ...  % VDE rhs is not vectorized
                      , 'stats'       , 'off'  );
   
end


methods
   
   % CONSTRUCTOR
   function obj = FPEX0_class_integration(varargin)
      if (nargin >= 1);   obj.integrator    = varargin{1};   end
      if (nargin >= 2);   obj.options       = varargin{2};   end
      if (nargin >= 3);   obj.VDEintegrator = varargin{3};   end
      if (nargin >= 4);   obj.VDEoptions    = varargin{4};   end
   end
   
   % Updater: Full Options
   function opts = updateOptions(obj, newOptions)
      obj.options = odeset(obj.options, newOptions);
      opts = obj.options;
   end
   
   % Updater: Jacobian only
   function opts = updateODEJacobian(obj, jacobianFcn)
      opts = obj.options;            % make copy of options
      opts.Jacobian = jacobianFcn;   % update JacobianFcn in copy
      obj.options = opts;            % save copy in options
   end

   
   
end

         
   
   
end