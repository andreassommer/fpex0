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
%
% INPUT:  integrator --> handle to Matlab integrator    (default: ode15s)
%         odeoptions --> integrator options from odeset (default: see source)
%
% OUTPUT: setup --> integator setup
%
% Author:  Andreas Sommer, Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%   
   

properties
   integrator = @ode15s
   options    = odeset( 'AbsTol'      , 1e-8     ...
                      , 'RelTol'      , 1e-6     ...
                      , 'BDF'         , 'off'    ...
                      , 'Vectorized'  , 'on'     ...
                      , 'stats'       , 'off'  );
end


methods
   
   % CONSTRUCTOR
   function obj = FPEX0_class_integrator(varargin)
      if (nargin >= 1);   obj.Integrator = varargin{1};   end
      if (nargin >= 2);   obj.Options    = varargin{2};   end
   end
   
   % Updater: Full Options
   function opts = updateOptions(obj, newOptions)
      obj.options = odeset(obj.options, newOptions);
      opts = obj.options;
   end
   
   % Updater: Jacobian only
   function opts = updateJacobian(obj, jacobianFcn)
      opts = obj.options;            % make copy of options
      opts.Jacobian = jacobianFcn;   % update JacobianFcn in copy
      obj.options = opts;            % save copy in options
   end

   
   
end

         
   
   
end