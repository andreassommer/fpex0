classdef FPEX0_class_grid < handle
   % class FPEX0_class_grid
   %
   % Class for storing the computation grid.
   %
   % CONSTRUCTOR:
   %
   %    FPEX0_class_grid(gridT, gridTdot)
   %
   %       gridT --> temperature grid (space grid for Fokker-Planck)
   %    gridTdot --> heating rate grid (time grid for Fokker-Planck)
   %
   %
   % Author:  Andreas Sommer, Sep2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %   
   
   properties (SetAccess = private)
      gridT       % x-grid = temperatures
      gridTdot    % t-grid = heating rates
   end
   
   properties (Dependent = true)
      N           % number of grid points
      h           % grid step size
   end
   

   methods
      
      % getter functions
      function N = get.N(obj)
         N = length(obj.gridT);
      end
      
      function h = get.h(obj)
         h = ( obj.gridT(end) - obj.gridT(1) ) / obj.N;
      end
      
      % CONSTRUCTOR
      function obj = FPEX0_class_grid(gridT, gridTdot)
         obj.gridT    = reshape( gridT   , [], 1);  % ensure it's a single column
         obj.gridTdot = reshape( gridTdot, [], 1);  % ensure it's a single column
      end
      
   end
   
   
end