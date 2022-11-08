classdef FPEX0_class_parameters < handle
   
   properties (SetAccess = private)
      p0             % initial values (first guess)
      p_lb           % lower bounds for parameters
      p_ub           % upper bounds for parameters
      idxFPdrift     % indices of parameter vector belonging to FP drift
      idxFPdiffusion % indices of parameter vector belonging to FP diffusion
      idxFPall       % indices of parameter vector belonging to all FP parameters
      idxIniDist     % indices of parameter vector belonging to initial distribution
   end
   
   properties (Dependent = true)
      count          % total number of parameters
   end
   

   methods

      %#ok<*MCSOH> % setter functions (silenced: not need to return obj in handle class)
      function obj = set.p0(obj,values);      obj.p0     = values;          end 
      function obj = set.p_lb(obj,values);    obj.p_lb   = values;          end 
      function obj = set.p_ub(obj,values);    obj.p_ub   = values;          end 
      
      % getter functions
      function val = get.count(obj);          val = length(obj.p0);         end

      % accessor functions for (external/input) parameter vector
      function val = extract_p_all(~,pp);           val = pp;                      end      
      function val = extract_p_FPdrift(obj,pp);     val = pp(obj.idxFPdrift);      end
      function val = extract_p_FPdiffusion(obj,pp); val = pp(obj.idxFPdiffusion);  end
      function val = extract_p_FPall(obj,pp);       val = pp(obj.idxFPall);        end
      function val = extract_p_iniDist(obj,pp);     val = pp(obj.idxIniDist);      end

      % CONSTRUCTOR
      function obj = FPEX0_class_parameters(  p0_FPdrift,   p0_FPdiffusion,   p0_iniDist, ...
                                            p_lb_FPdrift, p_lb_FPdiffusion, p_lb_iniDist, ...
                                            p_ub_FPdrift, p_ub_FPdiffusion, p_ub_iniDist)
         % quick error check on dimensions
         assert( isequal(size(p0_FPdrift)    , size(p_lb_FPdrift)    , size(p_ub_FPdrift))     );
         assert( isequal(size(p0_FPdiffusion), size(p_lb_FPdiffusion), size(p_ub_FPdiffusion)) );
         assert( isequal(size(p0_iniDist)    , size(p_lb_iniDist)    , size(p_ub_iniDist))     );
                                         
         % assemble the initial values, the lower bounds, the upper bounds
         obj.p0   = [ reshape(  p0_FPdrift,1,[]) , reshape(  p0_FPdiffusion,1,[]) , reshape(  p0_iniDist,1,[]) ];
         obj.p_lb = [ reshape(p_lb_FPdrift,1,[]) , reshape(p_lb_FPdiffusion,1,[]) , reshape(p_lb_iniDist,1,[]) ];
         obj.p_ub = [ reshape(p_ub_FPdrift,1,[]) , reshape(p_ub_FPdiffusion,1,[]) , reshape(p_ub_iniDist,1,[]) ];
                
         % setup the indices
         obj.idxFPdrift     = ( 1 : length(p0_FPdrift)     );
         obj.idxFPdiffusion = ( 1 : length(p0_FPdiffusion) ) + length(p0_FPdrift);
         obj.idxIniDist     = ( 1 : length(p0_iniDist)     ) + length(p0_FPdrift) + length(p0_FPdiffusion);
         obj.idxFPall       = [obj.idxFPdrift , obj.idxFPdiffusion];
                  
      end
      
      
   end
   
   
  

   
end