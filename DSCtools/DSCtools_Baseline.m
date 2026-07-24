classdef DSCtools_Baseline
% DSCTOOLS_BASELINE is a simple baseline description class containing function, onset, and offset.
%
% The function can be evaluated by calling the eval() function, or directly via its member
%
%
% Author:  Andreas Sommer, Jul2026
% code@andreas-sommer.eu
%

   properties
      blfun  {}
      onset  {mustBeNumeric}
      endset {mustBeNumeric}
      type   {}
      
   end
   methods
      function r = eval(obj, T)
         r = obj.blfun(T);
      end

      function obj = DSCtools_Baseline(blfun, onset, endset, type)
         obj.blfun  = blfun;
         obj.onset  = onset;
         obj.endset = endset;
         obj.type   = type;
      end

   end
end