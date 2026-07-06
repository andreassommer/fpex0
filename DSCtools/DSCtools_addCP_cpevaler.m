function val = DSCtools_addCP_cpevaler(T, cp, Tmin, Tmax)
% val = DSCtools_addCP_cpevaler(T, cp, Tmin, Tmax)
%
% Helper function to evaluate a cp function created by DSCtools_addCP().
%
%
% Author: Andreas Sommer, Jun2026
% code@andreas-sommer.eu

val = NaN(size(T));
nanIdx = (T<Tmin) | (T>Tmax);
val(~nanIdx) = ppval(cp.pp, T(~nanIdx));
if any(nanIdx)
   warning('Warning: leaving interpolation temperature range [%g, %g]', Tmin, Tmax)
end
  
end % of function
