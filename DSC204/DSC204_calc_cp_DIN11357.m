function cpS = DSC204_calc_cp_DIN11357(T, mS, dscS, cpR, mR, dscR, dsc0)
   % cpS = DSC204_calc_cp_DIN11357(T, mS, dscS, cpR, mR, dscR, dsc0) 
   %
   % Calculates the (apparent) specific heat capacity.
   %
   % It applies the "heat flow clibration" method: A known reference cp (of sapphire) 
   % is rescaled using the mass-ratio and signal ratio of reference and sample.
   % 
   %   cpS(T) = cpR(T) * mR/mS * (dscS(T)-dsc0(T)) / (dscR(T)-dsc0(T))
   %
   % INPUT:    T --> vector of temperatures to evaluate the cp-value  [double vector]
   %          mS --> mass of sample                                   [double]
   %        dscS --> DSC signal of sample (microvolts)                [function or double vector]
   %         cpR --> cp-values of reference                           [function or double vector]
   %          mR --> mass of reference                                [double]
   %        dscR --> DSC signal of reference (microvolts)             [function or double vector]
   %        dsc0 --> DSC signal with two empty crucibles              [function or double vector]
   %
   % OUTPUT: cpS --> cp-values of sample at specified temperatures.
   %
   % NOTES: * If vectors are given, they must coincide in size and meaning, 
   %          i.e.  dscS(k), dscR(k), dsc0(k), cpR(k) all correspond to temperature T(k).
   %        * If dscS and dscR are already baseline-corrected, set dsc0 to 0.
   %
   % DIN EN ISO 11357-4:2014
   % Plastics - Differential scanning calorimetry
   % Determination of specific heat capacity
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   % if some quantities are given as functions, generate the vectors
   if isa(dscS, 'function_handle'), dscS = dscS(T); end
   if isa(dscR, 'function_handle'), dscR = dscR(T); end
   if isa(dsc0, 'function_handle'), dsc0 = dsc0(T); end
   if isa(cpR , 'function_handle'), cpR  = cpR(T);  end
   
   % if dsc0 is missing, set it to zero
   if ~exist('dsc0', 'var'), dsc0 = 0; end
   
   % error check: dimensions
   assert(isequal(size(T), size(dscS), size(dscR), size(cpR)))
   assert(isequal(size(T), size(dsc0)) || isscalar(dsc0));
   
   % calculate cp
   cpS = cpR .* mR/mS .* (dscS-dsc0) ./ (dscR-dsc0);
   
end