classdef FPEX0_class_functions
   
   
   ERROR !! NOT TO BE USED !!
   
   properties
      
      FPdrift;     % function handle to FP drift
      FPdiffusion; % function handle to FP diffusion
      IniDist;     % function handle to Initial Distribution
   
   end
   
   
   methods
      
      % CONSTRUCTOR
      function obj = FPEX0_class_functions(driftFcn, diffusionFcn, iniDistFcn)
         obj.FPdrift     = driftFcn;
         obj.FPdiffusion = diffusionFcn;
         obj.IniDist     = iniDistFcn;
      end
      
   end
   
end