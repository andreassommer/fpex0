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
      N           % number of grid points
      h           % grid step size
   end
   

   methods
      
      % CONSTRUCTOR
      function obj = FPEX0_class_grid(gridT, gridTdot)
         obj.gridT    = reshape( gridT   , [], 1);  % ensure it's a single column
         obj.gridTdot = reshape( gridTdot, [], 1);  % ensure it's a single column
         obj.N        = length(obj.gridT);                            % number of grid points
         obj.h        = ( obj.gridT(end) - obj.gridT(1) ) / obj.N;    % grid step size
      end
      
   end
   
   
end